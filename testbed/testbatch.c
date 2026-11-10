/* General float batching regression. Distributed under the Little CMS MIT license. */
#include "lcms2_internal.h"
#include <fenv.h>
#include <stdio.h>
#include <stdint.h>

static unsigned Checks;
static char Context[128];
// Abort the regression test with the failed expression and its source line.
#define REQUIRE(x)                                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        if (!(x))                                                                                                      \
        {                                                                                                              \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #x);                                                       \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)
// Require bit-identical buffers and report the first differing byte with the active test context.
static void Equal(const void *a, const void *b, size_t bytes)
{
    if (memcmp(a, b, bytes))
    {
        size_t i;
        for (i = 0; i < bytes; ++i)
            if (((const unsigned char *)a)[i] != ((const unsigned char *)b)[i])
                break;
        fprintf(stderr, "%s comparison=%u byte=%zu actual=%02x expected=%02x\n", Context, Checks, i,
                ((const unsigned char *)a)[i], ((const unsigned char *)b)[i]);
        REQUIRE(0);
    }
    ++Checks;
}
// Float results must match exactly except for NaN encodings. Infinities and
// signed zero still require identical bits. Test padding uses finite sentinels.
static void EqualFloat(const void *a, const void *b, size_t bytes)
{
    size_t i;
    REQUIRE(bytes % sizeof(cmsFloat32Number) == 0);
    for (i = 0; i < bytes; i += sizeof(cmsFloat32Number)) {
        cmsUInt32Number x, y;
        memcpy(&x, (const unsigned char*)a + i, sizeof(x));
        memcpy(&y, (const unsigned char*)b + i, sizeof(y));
        if (x == y) continue;
        if ((x & 0x7fffffffu) > 0x7f800000u && (y & 0x7fffffffu) > 0x7f800000u) continue;
        fprintf(stderr, "%s float comparison=%u byte=%zu actual=%08x expected=%08x\n",
                Context, Checks, i, x, y);
        REQUIRE(0);
    }
    ++Checks;
}
// Generate reproducible finite, extended-range and exceptional float inputs for batch comparisons.
static float Value(unsigned i)
{
    static const cmsUInt32Number edge[] = {0,           0x80000000u, 1,           0x80000001u,
                                           0x7f800000u, 0xff800000u, 0x7fc12345u, 0x7f800001u};
    float f;
    if (i < sizeof(edge) / sizeof(edge[0]))
    {
        memcpy(&f, edge + i, sizeof(f));
        return f;
    }
    return ((int)(i * 37 % 401) - 50) / 300.0f;
}
// Compare a strided curve batch with scalar evaluation, including padding, then free the curve.
static void Curve(cmsToneCurve *c)
{
    float in[129 * 2], out[129 * 3], ref[129 * 3];
    unsigned i;
    REQUIRE(c && c->EvalBatch);
    for (i = 0; i < 129; ++i)
    {
        in[i * 2] = Value(i);
        in[i * 2 + 1] = 123;
    }
    for (i = 0; i < 129 * 3; ++i)
        out[i] = ref[i] = 321;
    c->EvalBatch(c, in, 2, out, 3, 0);
    Equal(out, ref, sizeof(out));
    for (i = 0; i < 129; ++i)
        ref[i * 3] = cmsEvalToneCurveFloat(c, in[i * 2]);
    c->EvalBatch(c, in, 2, out, 3, 129);
    Equal(out, ref, sizeof(out));
    cmsFreeToneCurve(c);
}
// Cover forward/inverse parametric, sampled and segmented curves under all supported rounding modes.
static void Curves(void)
{
    const int types[] = {1, 2, 3, 4, 5, 6, 7, 8, 108, 109};
    const int modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    unsigned t, m;
    int sign, saved = fegetround();
    for (m = 0; m < 4; ++m)
    {
        fesetround(modes[m]);
        for (t = 0; t < sizeof(types) / sizeof(types[0]); ++t)
            for (sign = -1; sign <= 1; sign += 2)
            {
                double p[] = {2.4, .9, .1, .08, .04, .01, .02, 0, 0, 0};
                snprintf(Context, sizeof(Context), "curve=%d mode=%u gamma=%g", sign * types[t], m, p[0]);
                Curve(cmsBuildParametricToneCurve(NULL, sign * types[t], p));
                p[0] = p[1] = p[3] = 0; // degenerate parameters retain original guards
                snprintf(Context, sizeof(Context), "curve=%d mode=%u gamma=%g", sign * types[t], m, p[0]);
                Curve(cmsBuildParametricToneCurve(NULL, sign * types[t], p));
            }
        {
            cmsUInt16Number table[] = {0, 1000, 35000, 65535};
            float samples[] = {0, .03f, .5f, 1};
            cmsCurveSegment segments[2] = {{0}};
            segments[0].x0 = -1e22f;
            segments[0].x1 = .5f;
            segments[0].Type = 1;
            segments[0].Params[0] = 1.8;
            segments[1].x0 = .5f;
            segments[1].x1 = 1e22f;
            segments[1].Type = 1;
            segments[1].Params[0] = 2.2;
            Curve(cmsBuildSegmentedToneCurve(NULL, 2, segments));
            Curve(cmsBuildTabulatedToneCurve16(NULL, 4, table));
            Curve(cmsBuildTabulatedToneCurveFloat(NULL, 4, samples));
        }
    }
    fesetround(saved);
}
// Compare scalar and batch stage outputs across block sizes and padding, then free the stage.
static void Stage(cmsStage *s, cmsBool Extended)
{
    const unsigned counts[] = {0, 1, 3, 4, 7, 8, 15, 16, 17, 63, 64, 65, 129};
    const unsigned stride = MAX_STAGE_CHANNELS;
    float *in = calloc(129 * stride, sizeof(float)), *out = calloc(129 * stride, sizeof(float)),
          *ref = calloc(129 * stride, sizeof(float));
    unsigned i, c, k;
    REQUIRE(s && s->EvalBatch && in && out && ref);
    for (i = 0; i < 129; ++i)
        for (c = 0; c < s->InputChannels; ++c)
            in[i * stride + c] = Extended ? Value(i + c) : ((i + c) % 97) / 96.0f;
    for (k = 0; k < sizeof(counts) / sizeof(counts[0]); ++k)
    {
        memset(out, 0xa5, 129 * stride * sizeof(float));
        memcpy(ref, out, 129 * stride * sizeof(float));
        for (i = 0; i < counts[k]; ++i)
            s->EvalPtr(in + i * stride, ref + i * stride, s);
        snprintf(Context, sizeof(Context), "stage=%x channels=%u/%u count=%u", s->Type, s->InputChannels,
                 s->OutputChannels, counts[k]);
        s->EvalBatch(in, stride, out, stride, counts[k], s);
        if (s->Type == cmsSigMatrixElemType) EqualFloat(out, ref, 129 * stride * sizeof(float));
        else Equal(out, ref, 129 * stride * sizeof(float));
    }
    free(in);
    free(out);
    free(ref);
    cmsStageFree(s);
}
// Create a deterministic test matrix with independent channel counts and optional offsets.
static cmsStage *Matrix(unsigned in, unsigned out, cmsBool offset)
{
    double *coeff = calloc(in * out, sizeof(double)), *off = calloc(out, sizeof(double));
    cmsStage *s;
    unsigned i;
    REQUIRE(coeff && off);
    for (i = 0; i < in * out; ++i)
        coeff[i] = (i % 7) / 11.0;
    for (i = 0; i < out; ++i)
        off[i] = (i % 5) / 13.0;
    s = cmsStageAllocMatrix(NULL, out, in, coeff, offset ? off : NULL);
    free(coeff);
    free(off);
    return s;
}
// Populate the test CLUT with a simple affine mapping for each of its three channels.
static int Sample(const float In[], float Out[], void *cargo)
{
    unsigned c;
    cmsUNUSED_PARAMETER(cargo);
    for (c = 0; c < 3; ++c)
        Out[c] = In[c] * .8f + .03f;
    return 1;
}
// Exercise matrix dimensions, offsets and the remaining built-in stages through the common batch checker.
static void Stages(void)
{
    unsigned i, j;
    const unsigned channels[] = {1, 3, 4, 7, 128};
    for (i = 0; i < 5; ++i)
        for (j = 0; j < 5; ++j)
        {
            Stage(Matrix(channels[i], channels[j], FALSE), TRUE);
            Stage(Matrix(channels[i], channels[j], TRUE), TRUE);
        }
    {
        double coeff[] = {NAN, 0, 0, 0, 1, 0, 0, 0, 1}, off[] = {0, 0, 0};
        Stage(cmsStageAllocMatrix(NULL, 3, 3, coeff, off), TRUE);
    }
    Stage(cmsStageAllocIdentity(NULL, 7), TRUE);
    Stage(_cmsStageAllocLab2XYZ(NULL), TRUE);
    Stage(_cmsStageAllocXYZ2Lab(NULL), TRUE);
    {
        cmsStage *s = cmsStageAllocCLutFloat(NULL, 3, 3, 3, NULL);
        REQUIRE(s);
        REQUIRE(cmsStageSampleCLutFloat(s, Sample, NULL, 0));
        Stage(s, FALSE);
        s = cmsStageAllocCLut16bit(NULL, 3, 3, 3, NULL);
        Stage(s, FALSE);
    }
}
// Removing the batch wrapper's finite scan must preserve the actual evaluators,
// including their exceptional inputs and CLUT bounds handling.
// Verify exceptional-input handling and CLUT bounds across rounding modes after removing wrapper scans.
static void StageInputChecks(void)
{
    const int modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    int saved = fegetround();
    cmsHPROFILE profile = cmsCreate_sRGBProfile();
    cmsHTRANSFORM h = cmsCreateTransform(profile, TYPE_RGB_FLT, profile, TYPE_RGB_FLT,
        INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | cmsFLAGS_NONEGATIVES);
    cmsStage* clip;
    unsigned m, dims, floating, i;
    REQUIRE(h);
    cmsCloseProfile(profile);
    for (clip = cmsPipelineGetPtrToFirstStage(cmsGetTransformPipeline(h)); clip; clip = cmsStageNext(clip))
        if (cmsStageType(clip) == cmsSigClipNegativesElemType) break;
    REQUIRE(clip);
    for (m = 0; m < sizeof(modes)/sizeof(modes[0]); ++m) {
        REQUIRE(fesetround(modes[m]) == 0);
        Stage(cmsStageAllocIdentity(NULL, 7), TRUE);
        Stage(_cmsStageAllocLab2XYZ(NULL), TRUE);
        Stage(_cmsStageAllocXYZ2Lab(NULL), TRUE);
        Stage(cmsStageDup(clip), TRUE);
        for (dims = 1; dims <= 4; ++dims)
            for (floating = 0; floating < 2; ++floating) {
                cmsStage* s = floating ? cmsStageAllocCLutFloat(NULL, 3, dims, 3, NULL) :
                                        cmsStageAllocCLut16bit(NULL, 3, dims, 3, NULL);
                _cmsStageCLutData* data;
                REQUIRE(s);
                data = (_cmsStageCLutData*)cmsStageData(s);
                for (i = 0; i < data->nEntries; ++i) {
                    if (floating) data->Tab.TFloat[i] = (i % 97) / 96.0f;
                    else data->Tab.T[i] = (cmsUInt16Number)((i * 751) % 65536);
                }
                Stage(s, TRUE);
            }
    }
    REQUIRE(fesetround(saved) == 0);
    cmsDeleteTransform(h);
}
static unsigned Events[8192], EventCount;
static cmsBool Recording;
// Append a callback event while recording is enabled so tests can check scalar callback ordering.
static void Record(unsigned event)
{
    if (Recording)
    {
        REQUIRE(EventCount < 8192);
        Events[EventCount++] = event;
    }
}
static _cmsStageEvalFn StockStage;
static cmsFormatterFloat StockIn, StockOut;
// Record a stage event before forwarding to the saved stock evaluator.
static void CustomStage(const float *in, float *out, const cmsStage *s)
{
    Record(2);
    StockStage(in, out, s);
}
// Record an input formatter event before invoking the saved stock formatter.
static cmsUInt8Number *CustomIn(_cmsTRANSFORM *p, float *values, cmsUInt8Number *buf, cmsUInt32Number stride)
{
    Record(1);
    return StockIn(p, values, buf, stride);
}
// Record an output formatter event before invoking the saved stock formatter.
static cmsUInt8Number *CustomOut(_cmsTRANSFORM *p, float *values, cmsUInt8Number *buf, cmsUInt32Number stride)
{
    Record(3);
    return StockOut(p, values, buf, stride);
}
// Record a parametric plugin call and return a simple formula for callback-order checks.
static double CustomCurve(int type, const double p[], double r)
{
    cmsUNUSED_PARAMETER(type);
    cmsUNUSED_PARAMETER(p);
    Record(4);
    return r * .7;
}
// Record a plugin interpolation call and fill each output channel deterministically.
static void CustomInterpolation(const float *in, float *out, const cmsInterpParams *p)
{
    unsigned c;
    Record(5);
    for (c = 0; c < p->nOutputs; ++c)
        out[c] = in[0] * .6f;
}
// Supply the recording float interpolator and leave unsupported interpolation formats unset.
static cmsInterpFunction InterpolationFactory(unsigned in, unsigned out, unsigned flags)
{
    cmsInterpFunction fn;
    cmsUNUSED_PARAMETER(in);
    cmsUNUSED_PARAMETER(out);
    fn.Lerp16 = NULL;
    if (flags & CMS_LERP_FLAGS_FLOAT)
        fn.LerpFloat = CustomInterpolation;
    return fn;
}
// Compare transform output with scalar stages and check strides, tails, in-place use and plugin ordering.
static void RunTransform(cmsHTRANSFORM handle, unsigned format, cmsBool plugin)
{
    _cmsTRANSFORM *p = (_cmsTRANSFORM *)handle;
    const cmsStage *stage;
    cmsBool hasMatrix = FALSE;
    const unsigned counts[] = {0, 1, 3, 4, 7, 8, 15, 16, 17, 63, 64, 65, 129};
    unsigned k, line, i, c;
    unsigned channels = T_CHANNELS(format) + T_EXTRA(format), bytes = T_BYTES(format) ? T_BYTES(format) : 8;
    unsigned row = 129 * channels * bytes + 64, plane = row * 3 + 64, size = plane * channels + 64;
    unsigned char *in = malloc(size), *out = malloc(size), *ref = malloc(size), *alias = malloc(size);
    unsigned expected[8192], n;
    REQUIRE(in && out && ref && alias);
    for (stage = p->Lut->Elements; stage != NULL; stage = stage->Next)
        if (stage->Type == cmsSigMatrixElemType) hasMatrix = TRUE;
    Recording = TRUE;
    for (i = 0; i < size; ++i)
        in[i] = 0;
    // All input formats in this test are floats; fill padded and planar regions too.
    for (i = 0; i + 4 <= size; i += 4)
    {
        float f = (i % 97) / 96.0f;
        memcpy(in + i, &f, 4);
    }
    for (k = 0; k < sizeof(counts) / sizeof(counts[0]); ++k)
    {
        snprintf(Context, sizeof(Context), "transform format=%x count=%u plugin=%d", format, counts[k], plugin);
        memset(ref, 0xa5, size);
        memcpy(out, ref, size);
        EventCount = 0;
        for (line = 0; line < 3; ++line)
        {
            cmsUInt8Number *from = in + line * row;
            cmsUInt8Number *to = ref + line * row;
            for (i = 0; i < counts[k]; ++i)
            {
                float a[cmsMAXCHANNELS] = {0}, b[cmsMAXCHANNELS] = {0};
                from = p->FromInputFloat(p, a, from, plane);
                cmsPipelineEvalFloat(a, b, p->Lut);
                to = p->ToOutputFloat(p, b, to, plane);
            }
        }
        n = EventCount;
        memcpy(expected, Events, n * sizeof(unsigned));
        EventCount = 0;
        // Alpha is handled outside the per-pixel callbacks.
        if (T_EXTRA(format))
            for (line = 0; line < 3; ++line)
                for (i = 0; i < counts[k]; ++i)
                    for (c = 3; c < channels; ++c)
                    {
                        size_t pos = line * row + (T_PLANAR(format) ? c * plane + i * 4 : (i * channels + c) * 4);
                        memcpy(ref + pos, in + pos, 4);
                    }
        cmsDoTransformLineStride(handle, in, out, counts[k], 3, row, row, plane, plane);
        if (hasMatrix) EqualFloat(ref, out, size); else Equal(ref, out, size);
        if (plugin)
        {
            REQUIRE(n == EventCount);
            Equal(expected, Events, n * sizeof(unsigned));
        }
        memcpy(alias, in, size);
        // Reference padding for the in-place call starts with the input buffer.
        for (line = 0; line < 3; ++line)
            for (i = 0; i < counts[k]; ++i)
                for (c = 0; c < channels; ++c)
                {
                    size_t pos = line * row + (T_PLANAR(format) ? c * plane + i * 4 : (i * channels + c) * 4);
                    memcpy(alias + pos, ref + pos, 4);
                }
        EventCount = 0;
        cmsDoTransformLineStride(handle, in, in, counts[k], 3, row, row, plane, plane);
        if (hasMatrix) EqualFloat(in, alias, size); else Equal(in, alias, size);
    }
    Recording = FALSE;
    free(in);
    free(out);
    free(ref);
    free(alias);
}
static cmsPipeline* ProposedPipeline;
// Transfer the proposed test pipeline into the transform during optimization, before batch preparation.
static cmsBool InstallPipeline(cmsPipeline** lut, cmsUInt32Number intent,
                               cmsUInt32Number* in, cmsUInt32Number* out, cmsUInt32Number* flags)
{
    cmsUNUSED_PARAMETER(intent);
    cmsUNUSED_PARAMETER(in);
    cmsUNUSED_PARAMETER(out);
    cmsUNUSED_PARAMETER(flags);
    if (ProposedPipeline == NULL) return FALSE;
    cmsPipelineFree(*lut);
    *lut = ProposedPipeline;
    ProposedPipeline = NULL;
    return TRUE;
}
// Create a test transform whose pipeline is supplied by the caller.
static cmsHTRANSFORM Transform(cmsContext ctx, unsigned format, cmsPipeline *lut)
{
    cmsPluginOptimization plugin = {{cmsPluginMagicNumber, LCMS_VERSION, cmsPluginOptimizationSig, NULL}, InstallPipeline};
    cmsHPROFILE profile = cmsCreate_sRGBProfileTHR(ctx);
    cmsHTRANSFORM h;
    // Install test stages during construction, before eligibility is prepared.
    REQUIRE(cmsPluginTHR(ctx, &plugin));
    ProposedPipeline = lut;
    h = cmsCreateTransformTHR(ctx, profile, format, profile, format, INTENT_RELATIVE_COLORIMETRIC,
                             cmsFLAGS_NOCACHE | cmsFLAGS_COPY_ALPHA);
    REQUIRE(h);
    REQUIRE(ProposedPipeline == NULL && ((_cmsTRANSFORM*)h)->Lut == lut);
    cmsCloseProfile(profile);
    return h;
}

// Compare each prepared formatter directly with its scalar reference, including
// untouched planes, padding and pointer advancement. The transform is created
// through the public API, so these also check that preparation actually occurs.
static void CheckFormatterBatch(_cmsTRANSFORM* p, cmsBool input)
{
    const unsigned counts[] = {0, 1, 15, 16, 17, 64};
    const int modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
    _cmsFloatFormatterBatch* layout = input ? &p->BatchInput : &p->BatchOutput;
    float values[64 * cmsMAXCHANNELS], reference[64 * cmsMAXCHANNELS];
    float buffer[80 * (cmsMAXCHANNELS + 2)], expected[80 * (cmsMAXCHANNELS + 2)];
    unsigned i, k, m;
    int saved = fegetround();
    REQUIRE(p->BatchLut == p->Lut && layout->Fn != NULL);
    for (m = 0; m < sizeof(modes)/sizeof(modes[0]); ++m) {
        REQUIRE(fesetround(modes[m]) == 0);
        for (k = 0; k < sizeof(counts)/sizeof(counts[0]); ++k) {
            cmsUInt8Number *end, *refend;
            for (i = 0; i < sizeof(values)/sizeof(values[0]); ++i)
                values[i] = reference[i] = input ? 321.0f : Value(i);
            for (i = 0; i < sizeof(buffer)/sizeof(buffer[0]); ++i)
                buffer[i] = expected[i] = input ? Value(i) : 123.0f;
            refend = (cmsUInt8Number*)expected;
            for (i = 0; i < counts[k]; ++i)
                refend = layout->Scalar(p, reference + i * cmsMAXCHANNELS, refend, 80 * sizeof(float));
            end = layout->Fn(p, values, (cmsUInt8Number*)buffer, 80 * sizeof(float), counts[k], cmsMAXCHANNELS);
            REQUIRE(end - (cmsUInt8Number*)buffer == refend - (cmsUInt8Number*)expected);
            snprintf(Context, sizeof(Context), "formatter input=%d format=%08x count=%u rounding=%d",
                     input, layout->Format, counts[k], modes[m]);
            Equal(values, reference, sizeof(values));
            Equal(buffer, expected, sizeof(buffer));
        }
    }
    REQUIRE(fesetround(saved) == 0);
}

// Check the actual plane layout, independently of the scalar/batch comparison.
// SWAPFIRST with no extra channels rotates the last color channel to the front.
static void PlanarSwapFirst(cmsHPROFILE link, unsigned format)
{
    const unsigned counts[] = {0, 1, 15, 16, 17, 64};
    enum { Rows = 3, Width = 64, Line = 71, Plane = Rows * Line + 7, Guard = 8 };
    float input[Rows * Width * cmsMAXCHANNELS];
    float output[Guard + cmsMAXCHANNELS * Plane + Guard];
    float expected[Guard + cmsMAXCHANNELS * Plane + Guard];
    const unsigned channels = T_CHANNELS(format);
    const float maximum = T_COLORSPACE(format) == PT_CMYK || T_COLORSPACE(format) >= PT_MCH5 ? 100.0f : 1.0f;
    unsigned scalar, swap, reverse, k, row, pixel, channel, i;

    for (row = 0; row < Rows; ++row)
        for (pixel = 0; pixel < Width; ++pixel)
            for (channel = 0; channel < channels; ++channel)
                input[(row * Width + pixel) * channels + channel] =
                    ((channel + 2 * pixel + row) % 9) * (maximum / 8.0f);

    for (scalar = 0; scalar < 2; ++scalar)
        for (swap = 0; swap < 2; ++swap)
            for (reverse = 0; reverse < 2; ++reverse) {
                unsigned out = format | PLANAR_SH(1) | SWAPFIRST_SH(1) | DOSWAP_SH(swap) | FLAVOR_SH(reverse);
                unsigned flags = scalar ? cmsFLAGS_NULLTRANSFORM : cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
                cmsHTRANSFORM h = cmsCreateTransform(scalar ? NULL : link, format, NULL, out,
                    INTENT_RELATIVE_COLORIMETRIC, flags);
                REQUIRE(h);
                if (!scalar) REQUIRE(((_cmsTRANSFORM*)h)->BatchLut == ((_cmsTRANSFORM*)h)->Lut);
                for (k = 0; k < sizeof(counts)/sizeof(counts[0]); ++k) {
                    for (i = 0; i < sizeof(output)/sizeof(output[0]); ++i)
                        output[i] = expected[i] = -99.0f;
                    for (row = 0; row < Rows; ++row)
                        for (pixel = 0; pixel < counts[k]; ++pixel)
                            for (channel = 0; channel < channels; ++channel) {
                                unsigned plane = swap ? channels - channel - 1 : channel;
                                float value = input[(row * Width + pixel) * channels + channel];
                                plane = (plane + 1) % channels;
                                expected[Guard + plane * Plane + row * Line + pixel] = reverse ? maximum - value : value;
                            }
                    cmsDoTransformLineStride(h, input, output + Guard, counts[k], Rows,
                        Width * channels * sizeof(float), Line * sizeof(float), 0, Plane * sizeof(float));
                    snprintf(Context, sizeof(Context), "planar SWAPFIRST scalar=%u format=%08x count=%u",
                        scalar, out, counts[k]);
                    Equal(output, expected, sizeof(output));
                }
                cmsDeleteTransform(h);
            }
}

// Cover stock formatter channel layouts and scalar adapters through transforms built by the public API.
static void FormatterBatches(void)
{
    const cmsColorSpaceSignature spaces[] = {cmsSigGrayData, cmsSigMCH2Data, cmsSigRgbData, cmsSigCmykData,
        cmsSigMCH5Data, cmsSigMCH6Data, cmsSigMCH7Data, cmsSigMCH8Data, cmsSigMCH9Data, cmsSigMCHAData,
        cmsSigMCHBData, cmsSigMCHCData, cmsSigMCHDData, cmsSigMCHEData, cmsSigMCHFData};
    const unsigned codes[] = {PT_GRAY, PT_MCH2, PT_RGB, PT_CMYK, PT_MCH5, PT_MCH6, PT_MCH7, PT_MCH8,
        PT_MCH9, PT_MCH10, PT_MCH11, PT_MCH12, PT_MCH13, PT_MCH14, PT_MCH15};
    cmsToneCurve* curves[cmsMAXCHANNELS];
    unsigned n, i, planar, swap, first, extra;
    curves[0] = cmsBuildGamma(NULL, 1.0);
    REQUIRE(curves[0]);
    for (i = 1; i < cmsMAXCHANNELS; ++i) curves[i] = curves[0];
    for (n = 1; n <= 15; ++n) {
        cmsHPROFILE link = cmsCreateLinearizationDeviceLink(spaces[n-1], curves);
        REQUIRE(link);
        PlanarSwapFirst(link, FLOAT_SH(1) | BYTES_SH(4) | CHANNELS_SH(n) | COLORSPACE_SH(codes[n-1]));
        for (planar = 0; planar < 2; ++planar)
            for (swap = 0; swap < 2; ++swap)
                for (first = 0; first < 2; ++first)
                    for (extra = 0; extra < 3; ++extra) {
                        unsigned format = FLOAT_SH(1) | BYTES_SH(4) | CHANNELS_SH(n) | COLORSPACE_SH(codes[n-1]) |
                            PLANAR_SH(planar) | DOSWAP_SH(swap) | SWAPFIRST_SH(first) | EXTRA_SH(extra);
                        cmsHTRANSFORM h = cmsCreateTransform(link, format, NULL, format | FLAVOR_SH(swap),
                            INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | cmsFLAGS_COPY_ALPHA);
                        REQUIRE(h);
                        CheckFormatterBatch((_cmsTRANSFORM*)h, TRUE);
                        CheckFormatterBatch((_cmsTRANSFORM*)h, FALSE);
                        cmsDeleteTransform(h);
                    }
        cmsCloseProfile(link);
    }
    cmsFreeToneCurve(curves[0]);
    // Lab normalization and premultiplied input (the latter uses the adapter).
    for (n = 0; n < 2; ++n)
        for (planar = 0; planar < 2; ++planar)
            for (extra = 0; extra < 3; ++extra) {
                cmsHPROFILE in = n ? cmsCreate_sRGBProfile() : cmsCreateLab4Profile(NULL);
                cmsHPROFILE out = cmsCreate_sRGBProfile();
                unsigned format = (n ? TYPE_RGB_FLT | PREMUL_SH(1) : TYPE_Lab_FLT) | PLANAR_SH(planar) | EXTRA_SH(extra);
                cmsHTRANSFORM h = cmsCreateTransform(in, format, out, TYPE_RGB_FLT,
                    INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
                REQUIRE(h);
                CheckFormatterBatch((_cmsTRANSFORM*)h, TRUE);
                cmsDeleteTransform(h);
                cmsCloseProfile(in); cmsCloseProfile(out);
            }
}
// Exceptional matrix arithmetic must preserve all non-NaN output bits.
static void MatrixExceptionalValues(void)
{
    const double values[] = {0.5, NAN, INFINITY, -INFINITY};
    unsigned i, offset;
    cmsStage *finite = Matrix(3, 3, TRUE);
    _cmsStageEvalBatchFn batch = finite->EvalBatch;
    cmsStageFree(finite);
    for (offset = 0; offset < 2; ++offset)
        for (i = 0; i < sizeof(values)/sizeof(values[0]); ++i) {
            double coeff[] = {1, 0, 0, 0, 1, 0, 0, 0, 1}, off[] = {0, 0, 0};
            cmsStage *matrix, *dup;
            cmsPipeline *lut = cmsPipelineAlloc(NULL, 3, 3), *copy;
            cmsHTRANSFORM h;
            if (offset) off[1] = values[i]; else coeff[0] = values[i];
            matrix = cmsStageAllocMatrix(NULL, 3, 3, coeff, off);
            REQUIRE(matrix && matrix->EvalBatch == batch);
            dup = cmsStageDup(matrix);
            REQUIRE(dup && dup->EvalBatch == matrix->EvalBatch);
            Stage(dup, TRUE);
            REQUIRE(lut && cmsPipelineInsertStage(lut, cmsAT_END, matrix));
            REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, Matrix(3, 3, FALSE)));
            copy = cmsPipelineDup(lut);
            REQUIRE(copy);
            h = Transform(NULL, TYPE_RGB_FLT, lut);
            RunTransform(h, TYPE_RGB_FLT, FALSE);
            cmsDeleteTransform(h);
            h = Transform(NULL, TYPE_RGB_FLT, copy);
            RunTransform(h, TYPE_RGB_FLT, FALSE);
            cmsDeleteTransform(h);
        }
    {
        cmsFloat64Number coeff[] = {INFINITY, -INFINITY, 0, 0, 1, 0, 0, 0, 1};
        cmsFloat64Number off[] = {0, 0, 0};
        const uint64_t payload1 = UINT64_C(0x7ff8000000001234);
        const uint64_t payload2 = UINT64_C(0xfff8000000005678);
        Stage(cmsStageAllocMatrix(NULL, 3, 3, coeff, off), TRUE);
        memcpy(&coeff[0], &payload1, sizeof(payload1));
        memcpy(&off[0], &payload2, sizeof(payload2));
        Stage(cmsStageAllocMatrix(NULL, 3, 3, coeff, off), TRUE);
    }
    // No prepared coefficient state: changed data and its duplicate use the same callback.
    {
        cmsStage *source = Matrix(3, 3, TRUE), *dup;
        ((_cmsStageMatrixData*)source->Data)->Double[0] = NAN;
        dup = cmsStageDup(source);
        REQUIRE(dup && dup->EvalBatch == source->EvalBatch);
        Stage(source, TRUE);
        Stage(dup, TRUE);
    }
    // A stage signature is not sufficient to identify matrix storage.
    {
        cmsStage *identity = cmsStageAllocIdentity(NULL, 3), *dup;
        identity->Type = cmsSigMatrixElemType;
        dup = cmsStageDup(identity);
        REQUIRE(dup && dup->EvalBatch == identity->EvalBatch);
        cmsStageFree(identity);
        Stage(dup, TRUE);
    }
}

// Exercise mixed built-in pipelines and verify that nested plugin callbacks keep scalar execution order.
static void Pipelines(void)
{
    unsigned mode, f;
    const unsigned formats[] = {TYPE_RGB_FLT, TYPE_RGBA_FLT, TYPE_RGB_FLT | PLANAR_SH(1), TYPE_RGBA_FLT | PLANAR_SH(1),
                                TYPE_BGR_FLT};
    for (mode = 0; mode < 4; ++mode)
        for (f = 0; f < 5; ++f)
        {
            cmsPipeline *lut = cmsPipelineAlloc(NULL, 3, 3);
            cmsHTRANSFORM h;
            if (mode == 1)
            {
                REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, Matrix(3, 128, TRUE)));
                REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, Matrix(128, 1, FALSE)));
                REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, Matrix(1, 3, TRUE)));
            }
            if (mode == 2)
            {
                cmsToneCurve *c[3];
                double p[] = {2.4, .9, .1, .08, .04, .01, .02};
                unsigned i;
                for (i = 0; i < 3; ++i)
                    c[i] = cmsBuildParametricToneCurve(NULL, i == 0 ? -5 : i == 1 ? 1 : 4, p);
                REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocToneCurves(NULL, 3, c)));
                cmsFreeToneCurveTriple(c);
            }
            if (mode == 3)
            {
                cmsStage *s = cmsStageAllocCLutFloat(NULL, 3, 3, 3, NULL);
                REQUIRE(cmsStageSampleCLutFloat(s, Sample, NULL, 0));
                REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, s));
            }
            {
                cmsPipeline *dup = cmsPipelineDup(lut);
                cmsPipelineFree(lut);
                lut = dup;
                REQUIRE(lut);
            }
            h = Transform(NULL, formats[f], lut);
            RunTransform(h, formats[f], FALSE);
            cmsDeleteTransform(h);
        }
    for (mode = 0; mode < 5; ++mode)
    {
        cmsContext ctx = cmsCreateContext(NULL, NULL);
        cmsPipeline *lut = cmsPipelineAlloc(ctx, 3, 3);
        cmsHTRANSFORM h;
        cmsPluginParametricCurves cp = {
            {cmsPluginMagicNumber, LCMS_VERSION, cmsPluginParametricCurveSig, NULL}, 1, {1}, {1}, CustomCurve};
        cmsPluginInterpolation ip = {{cmsPluginMagicNumber, LCMS_VERSION, cmsPluginInterpolationSig, NULL},
                                     InterpolationFactory};
        if (mode == 2)
        {
            cmsToneCurve *curves[3];
            unsigned c;
            for (c = 0; c < 3; ++c)
                curves[c] = cmsBuildGamma(ctx, 2.2);
            REQUIRE(cmsPluginTHR(ctx, &cp));
            REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocToneCurves(ctx, 3, curves)));
            cmsFreeToneCurveTriple(curves);
        }
        else if (mode == 4)
        {
            cmsToneCurve *curves[3];
            unsigned c;
            float samples[] = {0, .5f, 1};
            for (c = 0; c < 3; ++c)
                curves[c] = cmsBuildTabulatedToneCurveFloat(ctx, 3, samples);
            REQUIRE(cmsPluginTHR(ctx, &ip));
            REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocToneCurves(ctx, 3, curves)));
            cmsFreeToneCurveTriple(curves);
        }
        else if (mode == 3)
        {
            cmsStage *s = cmsStageAllocCLutFloat(ctx, 3, 3, 3, NULL), *dup;
            REQUIRE(s);
            REQUIRE(cmsPluginTHR(ctx, &ip));
            dup = cmsStageDup(s);
            cmsStageFree(s);
            REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, dup));
        }
        else
            REQUIRE(cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocIdentity(ctx, 3)));
        if (mode == 0)
        {
            StockStage = lut->Elements->EvalPtr;
            lut->Elements->EvalPtr = CustomStage;
        }
        h = Transform(ctx, TYPE_RGB_FLT, lut);
        if (mode == 1)
        {
            _cmsTRANSFORM *p = (_cmsTRANSFORM *)h;
            StockIn = p->FromInputFloat;
            StockOut = p->ToOutputFloat;
            p->FromInputFloat = CustomIn;
            p->ToOutputFloat = CustomOut;
        }
        RunTransform(h, TYPE_RGB_FLT, TRUE);
        cmsDeleteTransform(h);
        cmsDeleteContext(ctx);
    }
}

static unsigned Scheduled;
// Count scheduler invocations and forward the full request to the registered transform worker.
static void Scheduler(_cmsTRANSFORM* p, const void* in, void* out, cmsUInt32Number pixels,
                      cmsUInt32Number lines, const cmsStride* stride)
{
    _cmsTransform2Fn worker = _cmsGetTransformWorker(p);
    REQUIRE(worker);
    ++Scheduled;
    worker(p, in, out, pixels, lines, stride);
}

// Check prepared batching through a scheduler and scalar fallback after top-level state replacement.
static void BatchGuards(void)
{
    cmsPluginParalellization parallel = {{cmsPluginMagicNumber, LCMS_VERSION, cmsPluginParalellizationSig, NULL},
                                       2, 0, Scheduler};
    cmsContext ctx = cmsCreateContext(&parallel, NULL);
    cmsPipeline* lut = cmsPipelineAlloc(ctx, 3, 3);
    cmsPipeline* saved;
    cmsHTRANSFORM h;
    _cmsTRANSFORM* p;
    REQUIRE(ctx && lut && cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocIdentity(ctx, 3)));
    h = Transform(ctx, TYPE_RGB_FLT, lut);
    p = (_cmsTRANSFORM*)h;
    REQUIRE(p->BatchLut == p->Lut);
    RunTransform(h, TYPE_RGB_FLT, FALSE);
    REQUIRE(Scheduled != 0);

    // A replaced pipeline uses the scalar path; no execution-time cache writes.
    saved = p->Lut;
    p->Lut = cmsPipelineDup(saved);
    REQUIRE(p->Lut);
    StockStage = p->Lut->Elements->EvalPtr;
    p->Lut->Elements->EvalPtr = CustomStage;
    RunTransform(h, TYPE_RGB_FLT, TRUE);
    REQUIRE(p->BatchLut == saved);
    cmsDeleteTransform(h);
    cmsPipelineFree(saved);

    lut = cmsPipelineAlloc(ctx, 3, 3);
    REQUIRE(lut && cmsPipelineInsertStage(lut, cmsAT_END, cmsStageAllocIdentity(ctx, 3)));
    h = Transform(ctx, TYPE_RGB_FLT, lut);
    p = (_cmsTRANSFORM*)h;
    p->InputFormat = p->OutputFormat = TYPE_BGR_FLT;
    p->FromInputFloat = _cmsGetFormatter(ctx, TYPE_BGR_FLT, cmsFormatterInput, CMS_PACK_FLAGS_FLOAT).FmtFloat;
    p->ToOutputFloat = _cmsGetFormatter(ctx, TYPE_BGR_FLT, cmsFormatterOutput, CMS_PACK_FLAGS_FLOAT).FmtFloat;
    RunTransform(h, TYPE_BGR_FLT, FALSE);
    REQUIRE(p->BatchInput.Format == TYPE_RGB_FLT);
    cmsDeleteTransform(h);
    cmsDeleteContext(ctx);
}
// Run the curve, stage and transform batching regression groups and report the comparison count.
int main(void)
{
    Curves();
    Stages();
    StageInputChecks();
    Pipelines();
    MatrixExceptionalValues();
    FormatterBatches();
    BatchGuards();
    printf("PASS batch checks=%u\n", Checks);
    return 0;
}
