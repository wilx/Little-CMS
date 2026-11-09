/* General float batching regression. Distributed under the Little CMS MIT license. */
#include "lcms2_internal.h"
#include <fenv.h>
#include <stdio.h>

static unsigned Checks;
static char Context[128];
#define REQUIRE(x)                                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        if (!(x))                                                                                                      \
        {                                                                                                              \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #x);                                                       \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)
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
        Equal(out, ref, 129 * stride * sizeof(float));
    }
    free(in);
    free(out);
    free(ref);
    cmsStageFree(s);
}
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
static int Sample(const float In[], float Out[], void *cargo)
{
    unsigned c;
    cmsUNUSED_PARAMETER(cargo);
    for (c = 0; c < 3; ++c)
        Out[c] = In[c] * .8f + .03f;
    return 1;
}
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
        cmsStage *s = Matrix(3, 3, TRUE);
        ((_cmsStageMatrixData *)s->Data)->Double[0] = NAN;
        Stage(s, TRUE);
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
static unsigned Events[8192], EventCount;
static cmsBool Recording;
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
static void CustomStage(const float *in, float *out, const cmsStage *s)
{
    Record(2);
    StockStage(in, out, s);
}
static cmsUInt8Number *CustomIn(_cmsTRANSFORM *p, float *values, cmsUInt8Number *buf, cmsUInt32Number stride)
{
    Record(1);
    return StockIn(p, values, buf, stride);
}
static cmsUInt8Number *CustomOut(_cmsTRANSFORM *p, float *values, cmsUInt8Number *buf, cmsUInt32Number stride)
{
    Record(3);
    return StockOut(p, values, buf, stride);
}
static double CustomCurve(int type, const double p[], double r)
{
    cmsUNUSED_PARAMETER(type);
    cmsUNUSED_PARAMETER(p);
    Record(4);
    return r * .7;
}
static void CustomInterpolation(const float *in, float *out, const cmsInterpParams *p)
{
    unsigned c;
    Record(5);
    for (c = 0; c < p->nOutputs; ++c)
        out[c] = in[0] * .6f;
}
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
static void RunTransform(cmsHTRANSFORM handle, unsigned format, cmsBool plugin)
{
    _cmsTRANSFORM *p = (_cmsTRANSFORM *)handle;
    const unsigned counts[] = {0, 1, 3, 4, 7, 8, 15, 16, 17, 63, 64, 65, 129};
    unsigned k, line, i, c;
    unsigned channels = T_CHANNELS(format) + T_EXTRA(format), bytes = T_BYTES(format) ? T_BYTES(format) : 8;
    unsigned row = 129 * channels * bytes + 64, plane = row * 3 + 64, size = plane * channels + 64;
    unsigned char *in = malloc(size), *out = malloc(size), *ref = malloc(size), *alias = malloc(size);
    unsigned expected[8192], n;
    REQUIRE(in && out && ref && alias);
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
        Equal(ref, out, size);
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
        Equal(in, alias, size);
    }
    Recording = FALSE;
    free(in);
    free(out);
    free(ref);
    free(alias);
}
static cmsHTRANSFORM Transform(cmsContext ctx, unsigned format, cmsPipeline *lut)
{
    cmsHPROFILE profile = cmsCreate_sRGBProfileTHR(ctx);
    cmsHTRANSFORM h = cmsCreateTransformTHR(ctx, profile, format, profile, format, INTENT_RELATIVE_COLORIMETRIC,
                                            cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | cmsFLAGS_COPY_ALPHA);
    REQUIRE(h);
    cmsCloseProfile(profile);
    cmsPipelineFree(((_cmsTRANSFORM *)h)->Lut);
    ((_cmsTRANSFORM *)h)->Lut = lut;
    return h;
}
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
        h = Transform(ctx, TYPE_RGB_FLT, lut);
        if (mode == 0)
        {
            StockStage = lut->Elements->EvalPtr;
            lut->Elements->EvalPtr = CustomStage;
        }
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
int main(void)
{
    Curves();
    Stages();
    Pipelines();
    printf("PASS batch checks=%u\n", Checks);
    return 0;
}
