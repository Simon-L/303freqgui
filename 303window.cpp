#include "303window.hpp"

#include <cstring>
sst::surgext_rack::dsp::envelopes::ADAREnvelope vcf_env;
WowFilter wowFilter;

int reduced_size;
static constexpr int length = 10000;

float A = 2.14;
float B = 0.3;
float C = 0.449;
float D = 0.39;
float E = 4.602;

float fVco = 7.4;
float fVmod_amt = 0.887;
float fRes = 1.0;
float VaccMul = 2.547;

bool do_update = true;

float vcf_env_freq(float vcf_env, float Vco, float Vmod_amt, float Vacc, float VaccMul) {
    float Vmod_scale = 6.9*Vmod_amt+1.3;
    float Vmod_bias = -1.2*Vmod_amt+3;
    float Vmod = (Vmod_scale * vcf_env + Vmod_bias) - 3.2f; // 3.2 == Q9 bias
    // guest formula
    // Ic,11 = (A*Vco + B)*e^(C*Vmod + D*Vacc +E)
    return (A * Vco + B) * std::exp(C * Vmod + D * (Vacc * VaccMul) + E); // + D * Vacc
}

float plot_vcf_x[length];
float plot_vcf_y[length];
float plot_vcf_freq_maxvmod[length];
float plot_vcf_freq_minvmod[length];
float plot_vcf_freq_sweep[length];
float plot_vcf_vacc[length];

float plot_vcf_resvariation[20][length];

float curve_1_x[1000];
float curve_1_y[1000];

float maxFreqMinVmod;
float minFreqMinVmod;
float maxFreqMaxVmod;
float minFreqMaxVmod;

float fac0 = 0.125;
float fac1 = 0.2;
float fac2 = 0.25;
float fac3 = 0.325;
float fac4 = 0.875;
float bp0 = 0.0;
float bp1 = 0.25;
float bp2 = 0.5;
float bp3 = 0.75;
float bp4 = 1.0;

float xToFac(float x, float lo_bp, float hi_bp, float lo_fac, float hi_fac) {
    float Xmapping = 1.0f/(hi_bp - lo_bp);
    float mappedX = Xmapping * (x - lo_bp);
    return (hi_fac - lo_fac) * mappedX + lo_fac;
}

float xToY(float x) {
    if (x < bp1) return xToFac(x, bp0, bp1, fac0, fac1);
    if ((x >= bp1) && (x < bp2)) return xToFac(x, bp1, bp2, fac1, fac2);
    if ((x >= bp2) && (x < bp3)) return xToFac(x, bp2, bp3, fac2, fac3);
    if (x >= bp3) return xToFac(x, bp3, bp4, fac3, fac4);
}

void calculateCurve_1() {
    for (int i = 0; i < 1000; ++i)
    {
        float x = (float)i/1000.0;
        curve_1_x[i] = x;
        curve_1_y[i] = xToY(x);
    }
}

void calculateFreqRange() {
    vcf_env.immediatelyEnd();
    vcf_env.attackFrom(0.0f, 3, false, false);
    
    int reduced_x = 0;
    for (int samp = 0; samp < length; ++samp) {
        vcf_env.process(-9.482, -2.223, 3, 1, false);
        float env = vcf_env.output;
        if (samp % 35 == 0) {
            plot_vcf_x[reduced_x] = samp;
            plot_vcf_y[reduced_x] = env;
            float freq = vcf_env_freq(env, 7.4, 0.887, 0.0, 0.0);
            plot_vcf_freq_maxvmod[reduced_x] = freq;
            freq = vcf_env_freq(env, 7.4, 0.25, 0.0, 0.0);
            plot_vcf_freq_minvmod[reduced_x] = freq;
            reduced_x++;
        }
    }
    reduced_size = reduced_x;
    maxFreqMinVmod = vcf_env_freq(1.0016, 7.4, 0.25, 0.0, 0.0);
    minFreqMinVmod = vcf_env_freq(0.0, 7.4, 0.25, 0.0, 0.0);
    maxFreqMaxVmod = vcf_env_freq(1.0016, 7.4, 0.887, 0.0, 0.0);
    minFreqMaxVmod = vcf_env_freq(0.0, 7.4, 0.887, 0.0, 0.0);
    printf("Freq MinVmod: Max: %f Min: %f\n", maxFreqMinVmod, minFreqMinVmod);
    printf("Freq MaxVmod: Max: %f Min: %f\n", maxFreqMaxVmod, minFreqMaxVmod);
    // printf("reduced_size %d\n", reduced_size);
}

void calculateFreqSweep(bool accent) {
    vcf_env.immediatelyEnd();
    vcf_env.attackFrom(0.0f, 3, false, false);
    
    int reduced_x = 0;
    for (int samp = 0; samp < length; ++samp) {
        vcf_env.process(-9.482, -2.223, 3, 1, false);
        float Vacc = wowFilter.processSample(accent ? vcf_env.output * 9.91 : 0.0f);

        float env = vcf_env.output;
        if (samp % 35 == 0) {
            plot_vcf_x[reduced_x] = samp;
            float freq = vcf_env_freq(env, fVco, fVmod_amt, Vacc, VaccMul);
            freq = std::min(freq, 80000.0f);
            plot_vcf_freq_sweep[reduced_x] = freq;
            plot_vcf_vacc[reduced_x] = Vacc;
            reduced_x++;
        }
    }
    reduced_size = reduced_x;
    // printf("reduced_size %d\n", reduced_size);
}

void calculateFreqSweep2(bool accent, float* dest, float foo) {
    vcf_env.immediatelyEnd();
    vcf_env.attackFrom(0.0f, 3, false, false);
    
    int reduced_x = 0;
    for (int samp = 0; samp < length; ++samp) {
        vcf_env.process(-9.482, -2.223, 3, 1, false);
        float Vacc = wowFilter.processSample(accent ? vcf_env.output * (9.91 * foo) : 0.0f);

        float env = vcf_env.output;
        if (samp % 35 == 0) {
            float freq = vcf_env_freq(env, fVco, fVmod_amt, Vacc, VaccMul);
            // freq = std::min(freq, 80000.0f);
            dest[reduced_x] = freq;
            reduced_x++;
        }
    }
    reduced_size = reduced_x;
}

void init303Window() {
    vcf_env.activate(48000.0f);
    wowFilter.prepare(48000.0f);
    wowFilter.setResonancePot(fRes);
    calculateFreqRange();
    calculateCurve_1();
}

void show303Window() {
    ImGui::Begin("303 workbench", NULL, ImGuiWindowFlags_AlwaysAutoResize|ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoCollapse);

    if (ImGui::InputFloat("A - 2.14 - base scaling", &A, 0.01, 0.1)) {
        do_update |= true;
    }
    if (ImGui::InputFloat("B - 0.3 - base offset", &B, 0.01, 0.1)) {
        do_update |= true;
    }
    if (ImGui::InputFloat("C - 0.449 - exp vmod scaling", &C, 0.01, 0.1)) {
        do_update |= true;
    }
    if (ImGui::InputFloat("D - 0.39 - exp accent scaling", &D, 0.01, 0.1)) {
        do_update |= true;
    }
    if (ImGui::InputFloat("E - 4.602 - exponent offset", &E, 0.01, 0.1)) {
        do_update |= true;
    }

    if (ImGui::SliderFloat("Cutoff - 0.85 -> 7.4", &fVco, 0.5, 8.0)) {
        do_update |= true;

    }
    if (ImGui::SliderFloat("Envmod - 0.25 -> 0.887", &fVmod_amt, 0.2, 0.95)) {
        do_update |= true;
        
    }
    if (ImGui::SliderFloat("Resonance", &fRes, 0.0, 1.0)) {
        do_update |= true;
        wowFilter.setResonancePot(fRes);
    }
    if (ImGui::SliderFloat("VaccMul", &VaccMul, 0.0, 10.0)) {
        do_update |= true;
    }

    if (do_update) {
        calculateFreqRange();
        calculateFreqSweep(true);
        do_update = false;
    }


    if (ImGui::CollapsingHeader("Curve")) {
        ImGui::VSliderFloat("##fac0", ImVec2(15, 260), &fac0, 0.0, 1.0);
        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            ImGui::SetTooltip("fac0 %.3f", fac0);
        }
        ImGui::SameLine();
        ImGui::VSliderFloat("##fac1", ImVec2(15, 260), &fac1, 0.0, 1.0);
        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            ImGui::SetTooltip("fac1 %.3f", fac1);
        }
        ImGui::SameLine();
        ImGui::VSliderFloat("##fac2", ImVec2(15, 260), &fac2, 0.0, 1.0);
        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            ImGui::SetTooltip("fac2 %.3f", fac2);
        }
        ImGui::SameLine();
        ImGui::VSliderFloat("##fac3", ImVec2(15, 260), &fac3, 0.0, 1.0);
        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            ImGui::SetTooltip("fac3 %.3f", fac3);
        }
        ImGui::SameLine();
        ImGui::VSliderFloat("##fac4", ImVec2(15, 260), &fac4, 0.0, 1.0);
        if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
            ImGui::SetTooltip("fac4 %.3f", fac4);
        }
        calculateCurve_1();
        ImGui::SameLine();
        if (ImPlot::BeginPlot("Curve 1")) {
            ImPlot::SetupLegend(ImPlotLocation_North|ImPlotLocation_East, 0);
            ImPlot::SetupAxis(ImAxis_Y1, "", ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("Curve 1", curve_1_x, curve_1_y, 1000);
            ImPlot::EndPlot();
        }
    }

    wowFilter.prepare(48000);
    if (ImPlot::BeginPlot("Res variation", ImVec2(-1.0,-1.0))) {
        ImPlot::SetupLegend(ImPlotLocation_North|ImPlotLocation_East, 0);
        ImPlot::SetupAxis(ImAxis_Y1, "", ImPlotAxisFlags_AutoFit);
        for (int i = 0; i < 20; ++i)
        {   
            char plot_name[16];
            std::sprintf(plot_name, "Res %.3f", 1.0 - (0.05 * i));
            wowFilter.setResonancePot(1.0 - (0.05 * i));
            calculateFreqSweep2(true, plot_vcf_resvariation[i], xToY(1.0 - (0.05 * i)));
            ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4((0.05 * i),(0.05 * i),(0.05 * i),1.0f));
            ImPlot::PlotLine(plot_name, plot_vcf_x, plot_vcf_resvariation[i], length);
            ImPlot::PopStyleColor();
        }
        ImPlot::EndPlot();
    }
    wowFilter.prepare(48000);
    wowFilter.setResonancePot(fRes);

    if (ImGui::CollapsingHeader("Frequency ranges plot")) { 
        ImGui::Text("Max freq for max Vmod: %.2fHz", maxFreqMaxVmod); ImGui::SameLine();
        ImGui::TextColored(ImVec4(0.1372549086809158f, 0.1568627506494522f, 0.1098039224743843f, 1.0f), "|"); ImGui::SameLine();
        ImGui::Text("Min freq for max Vmod: %.2fHz", minFreqMaxVmod);
        ImGui::Separator();
        ImGui::Text("Max freq for min Vmod: %.2fHz", maxFreqMinVmod); ImGui::SameLine();
        ImGui::TextColored(ImVec4(0.1372549086809158f, 0.1568627506494522f, 0.1098039224743843f, 1.0f), "|"); ImGui::SameLine();
        ImGui::Text("Min freq for min Vmod: %.2fHz", minFreqMinVmod);
        if (ImPlot::BeginPlot("Frequency range")) {
            ImPlot::SetupLegend(ImPlotLocation_North|ImPlotLocation_East, 0);
            ImPlot::SetupAxis(ImAxis_Y1, "Frequency (Hz)", ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("Freq min vmod", plot_vcf_x, plot_vcf_freq_minvmod, reduced_size);
            ImPlot::PlotLine("Freq max vmod", plot_vcf_x, plot_vcf_freq_maxvmod, reduced_size);
            ImPlot::EndPlot();
        }
    }

    if (ImGui::CollapsingHeader("Frequency sweep plot")) { 
        if (ImPlot::BeginPlot("Frequency sweep")) {
            ImPlot::SetupLegend(ImPlotLocation_North|ImPlotLocation_East, 0);
            ImPlot::SetupAxis(ImAxis_Y1, "Frequency (Hz)", ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxis(ImAxis_Y2, "Voltage (V)", ImPlotAxisFlags_AutoFit);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
            ImPlot::PlotLine("Filter frequency sweep", plot_vcf_x, plot_vcf_freq_sweep, reduced_size);
            ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
            ImPlot::PlotLine("Vacc", plot_vcf_x, plot_vcf_vacc, reduced_size);
            ImPlot::EndPlot();
        }
    }

    ImGui::End();
}