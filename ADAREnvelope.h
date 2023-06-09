// Taken from Surge XT for Rack, minimal changes - 2023 SL
// Original license below:
/*
 * SurgeXT for VCV Rack - a Surge Synth Team product
 *
 * Copyright 2019 - 2022, Various authors, as described in the github
 * transaction log.
 *
 * SurgeXT for VCV Rack is released under the Gnu General Public Licence
 * V3 or later (GPL-3.0-or-later). The license is found in the file
 * "LICENSE" in the root of this repository or at
 * https://www.gnu.org/licenses/gpl-3.0.en.html
 *
 * All source for Surge XT for VCV Rack is available at
 * https://github.com/surge-synthesizer/surge-rack/
 */

#include <cmath>
#include <algorithm>

#ifndef SXT_RACK_ADARENVELOPE_H
#define SXT_RACK_ADARENVELOPE_H

namespace sst::surgext_rack::dsp::envelopes
{
struct ADAREnvelope
{
    static constexpr int tuning_table_size = 512;
    float table_envrate_linear alignas(16)[512];
    double dsamplerate_os{0};

    // Constants from Surge
    static const int BLOCK_SIZE = 2;
    const int BLOCK_SIZE_OS = BLOCK_SIZE * 2; // default oversampling is 2
    const float BLOCK_SIZE_INV = (1.f / BLOCK_SIZE);

    float sample_rate = 0.0;

    ADAREnvelope()
    {   
        for (int i = 0; i < BLOCK_SIZE; ++i) {
            outputCache[i] = 0;
        }
    }

    void activate(float sr)
    {
        sample_rate = sr;
        dsamplerate_os = sample_rate * 2.0;
        // from SurgeStorage::init_tables
        for (int i = 0; i < tuning_table_size; i++)
        {
            double k = dsamplerate_os * pow(2.0, (((double)i - 256.0) / 16.0)) / (double)BLOCK_SIZE_OS;
            table_envrate_linear[i] = (float)(1.f / k);
        }
    }

    // from SurgeStorage
    float envelope_rate_linear(float x)
    {
        x *= 16.f;
        x += 256.f;
        int e = (int)x;
        float a = x - (float)e;

        return (1 - a) * table_envrate_linear[e & 0x1ff] + a * table_envrate_linear[(e + 1) & 0x1ff];
    }

    enum Stage
    {
        s_attack,
        s_hold,
        s_decay,
        s_analog_residual_decay,
        s_eoc,
        s_complete
    } stage{s_complete};

    bool isDigital{true};
    bool isGated{false};

    float phase{0}, start{0};

    float output{0}, eoc_output{0};
    float outputCache[BLOCK_SIZE], outBlock0{0.f};
    int current{BLOCK_SIZE};
    int eoc_countdown{0};

    // Analog Mode
    float v_c1{0}, v_c1_delayed{0.f};
    bool discharge{false};

    void attackFrom(float fv, int ashp, bool id, bool ig)
    {
        float f = fv;

        switch (ashp)
        {
        case 0:
            // target = sqrt(target);
            f = f * f;
            break;
        case 2:
            // target = target * target * target;
            f = pow(f, 1.0 / 3.0);
            break;
        case 3:
            f = f; // ignored
            break;
        }
        phase = f;
        stage = s_attack;
        current = BLOCK_SIZE;
        outBlock0 = f;
        isDigital = id;
        isGated = ig;
        eoc_output = 0;
        eoc_countdown = 0;

        v_c1 = f;
        v_c1_delayed = f;
        discharge = false;
    }

    template <bool gated>
    inline float stepDigital(const float a, const float d, const bool gateActive)
    {
        float target = 0;
        switch (stage)
        {
        default:
            break;
        case s_attack:
        {
            phase += envelope_rate_linear(a);
            if (phase >= 1)
            {
                phase = 1;
                if constexpr (gated)
                    stage = s_hold;
                else
                    stage = s_decay;
            }
            if constexpr (gated)
                if (!gateActive)
                    stage = s_decay;
            target = phase;
            break;
        }
        case s_decay:
        {
            phase -= envelope_rate_linear(d);
            if (phase <= 0)
            {
                phase = 0;
                stage = s_eoc;
                eoc_countdown = (int)std::round(sample_rate * 0.01);
            }
            target = phase;
        }
        }
        return target;
    }

    inline void immediatelyEnd()
    {
        stage = s_complete;
        eoc_output = 0;
        output = 0;
        eoc_countdown = 0;
        current = BLOCK_SIZE;
        outBlock0 = 0;
    }
    inline void process(const float a, const float d, const int ashape, const int dshape,
                        const bool gateActive)
    {
        if (stage == s_complete)
        {
            output = 0;
            return;
        }

        if (stage == s_eoc)
        {
            output = 0;
            eoc_output = 1;

            eoc_countdown--;
            if (eoc_countdown == 0)
            {
                eoc_output = 0;
                stage = s_complete;
            }
            return;
        }
        eoc_output = 0;

        if (stage == s_analog_residual_decay && eoc_countdown)
        {
            eoc_output = 1;
            eoc_countdown--;
        }
        if (current == BLOCK_SIZE)
        {
            float target = 0;
            if (isGated && stage == s_hold)
            {
                target = 1;

                if (!gateActive)
                {
                    phase = 1;
                    stage = s_decay;
                }
            }
            else if (isDigital)
            {
                if (isGated)
                    target = stepDigital<true>(a, d, gateActive);
                else
                    target = stepDigital<false>(a, d, gateActive);
            }
            else
            {
                const float coeff_offset = 2.f - std::log2(sample_rate * BLOCK_SIZE_INV);

                auto ndc = (v_c1_delayed >= 0.99999f);
                if (ndc && !discharge)
                {
                    phase = 1;
                    stage = isGated ? s_hold : s_decay;
                }

                if (isGated && !discharge)
                {
                    ndc = !gateActive;
                }
                discharge = ndc || discharge;
                v_c1_delayed = v_c1;

                static constexpr float v_gate = 1.02f;
                auto v_attack = (!discharge) * v_gate;
                auto v_decay = (!discharge) * v_gate;

                // In this case we only need the coefs in their stage
                float coef_A = !discharge ? powf(2.f, std::min(0.f, coeff_offset - a)) : 0;
                float coef_D = discharge ? powf(2.f, std::min(0.f, coeff_offset - d)) : 0;

                auto diff_v_a = std::max(0.f, v_attack - v_c1);
                auto diff_v_d = std::min(0.f, v_decay - v_c1);

                v_c1 = v_c1 + diff_v_a * coef_A + diff_v_d * coef_D;
                target = v_c1;

                if (stage == s_decay)
                {
                    phase -= envelope_rate_linear(d);
                    if (phase <= 0)
                    {
                        eoc_countdown = (int)std::round(sample_rate * 0.01);
                        stage = s_analog_residual_decay;
                    }
                }
                if (v_c1 < 1e-6 && discharge)
                {
                    v_c1 = 0;
                    v_c1_delayed = 0;
                    discharge = false;
                    target = 0;
                    if (stage != s_analog_residual_decay)
                    {
                        eoc_countdown = (int)std::round(sample_rate * 0.01);
                        stage = s_eoc;
                    }
                    else
                    {
                        stage = s_complete;
                        eoc_countdown = 0;
                    }
                }
            }

            if (stage == s_attack)
            {
                switch (ashape)
                {
                case 0:
                    target = sqrt(target);
                    break;
                case 2:
                    target = target * target * target;
                    break;
                case 3:
                    target = pow(sqrt(target), 0.35);
                    break;
                }
            }
            else
            {
                switch (dshape)
                {
                case 0:
                    target = sqrt(target);
                    break;
                case 2:
                    target = target * target * target;
                    break;
                }
            }

            float dO = (target - outBlock0) * BLOCK_SIZE_INV;
            for (int i = 0; i < BLOCK_SIZE; ++i)
            {
                outputCache[i] = outBlock0 + dO * i;
            }
            outBlock0 = target;

            current = 0;
        }

        output = outputCache[current];
        current++;
    }
};
} // namespace sst::surgext_rack::dsp::envelopes
#endif // RACK_HACK_ADARENVELOPE_H
