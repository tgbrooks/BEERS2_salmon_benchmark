sample_ids = [1,2]
lanes_used = [1]

run_configs = {
    'unbiased': {
        'GC_bias': 'none',
        'pos_3prime_bias': 'none',
        'primer_bias': 'none',
    },
    'GC_bias_med': {
        'GC_bias': 'med',
        'pos_3prime_bias': 'none',
        'primer_bias': 'none',
    },
    'GC_bias_high': {
        'GC_bias': 'high',
        'pos_3prime_bias': 'none',
        'primer_bias': 'none',
    },
    'pos_3prime_bias_med': {
        'GC_bias': 'none',
        'pos_3prime_bias': 'med',
        'primer_bias': 'none',
    },
    'pos_3prime_bias_high': {
        'GC_bias': 'none',
        'pos_3prime_bias': 'high',
        'primer_bias': 'none',
    },
    'primer_bias_high': {
        'GC_bias': 'none',
        'pos_3prime_bias': 'none',
        'primer_bias': 'high',
    },
    'all_bias': {
        'GC_bias': 'high',
        'pos_3prime_bias': 'high',
        'primer_bias': 'high',
    },
}

GC_bias = dict(
    none = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": 0.0,
    },
    med = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -10,
    },
    high = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -50,
    },
)

pos_3prime_bias = dict(
    none = {
        "breakpoint_prob_per_base": 0.0,
    },
    med = {
        "breakpoint_prob_per_base": 0.0002,
    },
    high = {
        "breakpoint_prob_per_base": 0.001,
    },
)

primer_bias = dict(
    none = {
        "perfect_priming": "true",
    },
    high = {
        "perfect_priming": "false",
    },
)

# Insert the BEERS config details into each of run configs
# according to the desired bias values
for run, config in run_configs.items():
    config.update(GC_bias[config['GC_bias']])
    config.update(pos_3prime_bias[config['pos_3prime_bias']])
    config.update(primer_bias[config['primer_bias']])
