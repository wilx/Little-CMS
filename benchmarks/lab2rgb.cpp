// LCMS-only benchmark using captured Lab rows and randomized library order.
#include "clock.hpp"
#include <lcms2.h>
#include <lcms2_plugin.h>
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _WIN32
#ifndef NOMINMAX
// Do not import Windows min/max macros into the benchmark.
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <memory>

// Resolve an LCMS entry point from the selected library rather than link-time defaults.
template<class T> T symbol(void* library, const char* name) {
#ifdef _WIN32
    const auto address = GetProcAddress(static_cast<HMODULE>(library), name);
    if (!address) { std::fprintf(stderr, "%s: Windows error %lu\n", name, GetLastError()); std::exit(2); }
    // Win32 returns a generic function pointer. Preserve its representation
    // without a compiler warning about calling through an incompatible type.
    T fn;
    static_assert(sizeof(fn) == sizeof(address), "Windows function pointer size mismatch");
    std::memcpy(&fn, &address, sizeof(fn));
#else
    auto fn = reinterpret_cast<T>(dlsym(library, name));
    if (!fn) { std::fprintf(stderr, "%s\n", dlerror()); std::exit(2); }
#endif
    return fn;
}
struct Library {
    void* handle;
// Declare each public LCMS entry point with the type from the LCMS headers.
#define API(name) decltype(&::name) name;
    API(cmsCreateLab4Profile) API(cmsOpenProfileFromFile) API(cmsCloseProfile)
    API(cmsCreateTransform) API(cmsDeleteTransform) API(cmsDoTransform)
    API(cmsGetEncodedCMMversion) API(cmsGetTransformPipeline)
    API(cmsPipelineGetPtrToFirstStage) API(cmsStageNext) API(cmsStageType)
    API(cmsStageInputChannels) API(cmsStageOutputChannels) API(cmsStageData)
#undef API
    // Keep each measured library loaded until its transforms have been destroyed.
    explicit Library(const char* path) {
#ifdef _WIN32
        handle = LoadLibraryA(path);
        if (!handle) { std::fprintf(stderr, "%s: Windows error %lu\n", path, GetLastError()); std::exit(2); }
#else
        handle = dlopen(path, RTLD_NOW | RTLD_LOCAL);
        if (!handle) { std::fprintf(stderr, "%s\n", dlerror()); std::exit(2); }
#endif
#define API(name) name = symbol<decltype(name)>(handle, #name);
        API(cmsCreateLab4Profile) API(cmsOpenProfileFromFile) API(cmsCloseProfile)
        API(cmsCreateTransform) API(cmsDeleteTransform) API(cmsDoTransform)
        API(cmsGetEncodedCMMversion) API(cmsGetTransformPipeline)
        API(cmsPipelineGetPtrToFirstStage) API(cmsStageNext) API(cmsStageType)
        API(cmsStageInputChannels) API(cmsStageOutputChannels) API(cmsStageData)
#undef API
        if (cmsGetEncodedCMMversion() != 2200) std::exit(2);
    }
    // Release the native module handle after all calls through it have finished.
    ~Library() {
#ifdef _WIN32
        FreeLibrary(static_cast<HMODULE>(handle));
#else
        dlclose(handle);
#endif
    }
};
constexpr size_t width=7370, rows=32, pixels=width*rows;
// Read one complete captured fixture, rejecting truncated or unexpected layouts.
static std::vector<float> load(const std::string& path) {
    std::ifstream file(path,std::ios::binary|std::ios::ate);
    if (!file || file.tellg()!=std::streampos(pixels*3*sizeof(float))) std::exit(2);
    std::vector<float> data(pixels*3);
    file.seekg(0);file.read(reinterpret_cast<char*>(data.data()),data.size()*sizeof(float));
    if (!file) std::exit(2);
    return data;
}
// Count differing float bit patterns, independently of numerical tolerances.
static size_t differences(const std::vector<float>& a,const std::vector<float>& b) {
    size_t count=0;
    for(size_t i=0;i<a.size();++i) if(std::memcmp(&a[i],&b[i],sizeof(float))) ++count;
    return count;
}
// Time transform execution; input loading, construction and validation are outside timers.
int main(int argc,char** argv) {
    if(argc!=5 && argc!=6) {std::fprintf(stderr,"usage: %s libraries.tsv fixtures rounds seed [passes=1]\n",argv[0]);return 2;}
    std::ifstream list(argv[1]);
    std::vector<std::string> names;
    std::vector<std::unique_ptr<Library>> libraries;
    std::string line;
    while(std::getline(list,line)) {
        if(!line.empty() && line.back()=='\r')line.pop_back();
        if(line.empty())continue;
        const auto tab=line.find('\t');
        if(tab==std::string::npos || tab==0 || tab+1==line.size())return 2;
        names.push_back(line.substr(0,tab));
        libraries.emplace_back(new Library(line.substr(tab+1).c_str()));
    }
    if(libraries.empty())return 2;
    int repeats=std::atoi(argv[3]);if(repeats<1)return 2;
    const int passes=argc==6 ? std::atoi(argv[5]) : 1;
    if(passes<1)return 2;
    const double measured_pixels=double(pixels)*passes;
    auto input=load(std::string(argv[2])+"/input-lab-f32.bin");
    auto reference=load(std::string(argv[2])+"/expected-rgb-f32.bin");
    std::vector<float> output(input.size());
    std::vector<cmsHTRANSFORM> transforms(libraries.size());
    for(size_t v=0;v<libraries.size();++v) {
        auto& l=*libraries[v];
        auto in=l.cmsCreateLab4Profile(nullptr);
        auto out=l.cmsOpenProfileFromFile((std::string(argv[2])+"/RTv4_sRGB.icc").c_str(),"r");
        if(!in || !out)return 2;
        transforms[v]=l.cmsCreateTransform(in,TYPE_Lab_FLT,out,TYPE_RGB_FLT,INTENT_RELATIVE_COLORIMETRIC,
            cmsFLAGS_NOOPTIMIZE|cmsFLAGS_NOCACHE|cmsFLAGS_BLACKPOINTCOMPENSATION);
        l.cmsCloseProfile(in);l.cmsCloseProfile(out);
        if(!transforms[v])return 2;
        auto pipeline=l.cmsGetTransformPipeline(transforms[v]);if(!pipeline)return 2;
        auto stage=l.cmsPipelineGetPtrToFirstStage(pipeline);
        const cmsStageSignature expected[]={cmsSigLab2XYZElemType,cmsSigMatrixElemType,cmsSigCurveSetElemType};
        for(int i=0;i<3;++i) {
            if(!stage || l.cmsStageType(stage)!=expected[i] || l.cmsStageInputChannels(stage)!=3 || l.cmsStageOutputChannels(stage)!=3)return 2;
            if(i==1 && static_cast<_cmsStageMatrixData*>(l.cmsStageData(stage))->Offset)return 2;
            stage=l.cmsStageNext(stage);
        }
        if(stage)return 2;
        std::fprintf(stderr,"PIPELINE %s: Lab2XYZ -> Matrix3x3(no offset) -> Curves\n",names[v].c_str());
    }
    omp_set_dynamic(0);
    std::mt19937 random(unsigned(std::strtoul(argv[4],nullptr,10)));
    std::puts("round,threads,revision,seconds,ns_per_pixel,process_cpu_seconds,worker_cpu_seconds,worker_cpu_ns_per_pixel,passes");
    for(int round=-3;round<repeats;++round) {
        int counts[]={1,4};std::shuffle(counts,counts+2,random);
        for(int threads:counts) {
            omp_set_num_threads(threads);
            std::vector<size_t> order(libraries.size());
            for(size_t v=0;v<order.size();++v)order[v]=v;
            std::shuffle(order.begin(),order.end(),random);
            for(size_t v:order) {
                std::fill(output.begin(),output.end(),-12345.0f);
                double worker_cpu=0;
                const auto begin=bench::WallClock::now();
                const auto process_begin=bench::ProcessCpuClock::now();
                #pragma omp parallel reduction(+:worker_cpu)
                {
                    const auto thread_begin=bench::ThreadCpuClock::now();
                    for(int pass=0;pass<passes;++pass) {
                        // Identical static schedules in this parallel region assign
                        // the same rows to each worker, so passes need no barrier.
                        #pragma omp for schedule(static) nowait
                        for(int row=0;row<int(rows);++row)
                            libraries[v]->cmsDoTransform(transforms[v],input.data()+row*width*3,output.data()+row*width*3,width);
                    }
                    worker_cpu+=bench::seconds(bench::ThreadCpuClock::now()-thread_begin);
                }
                const double process_cpu=bench::seconds(bench::ProcessCpuClock::now()-process_begin);
                const double seconds=bench::seconds(bench::WallClock::now()-begin);
                if(std::memcmp(reference.data(),output.data(),output.size()*sizeof(float))) {
                    std::fprintf(stderr,"FAIL %s threads=%d float differences=%zu\n",names[v].c_str(),threads,differences(reference,output));return 1;
                }
                if(round>=0)std::printf("%d,%d,%s,%.9f,%.6f,%.9f,%.9f,%.6f,%d\n",round,threads,names[v].c_str(),seconds,seconds*1e9/measured_pixels,process_cpu,worker_cpu,worker_cpu*1e9/measured_pixels,passes);
            }
        }
    }
    for(size_t v=0;v<libraries.size();++v)libraries[v]->cmsDeleteTransform(transforms[v]);
    std::fprintf(stderr,"PASS %zu variants, %d rounds, 1/4 threads, every output bit exact\n",libraries.size(),repeats);
}
