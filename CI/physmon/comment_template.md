## :bar_chart: Physics performance monitoring for {{ commit }}
{% if has_errors %}
> :red_square: **ERROR** The result has missing elements! 
> This is likely a physmon job failure
{% endif %}

[Full report]({{ url }}/)
CKF: {{ make_url("seeded", "ckf_seeded.html") }}, {{ make_url("truth smeared", "ckf_truth_smeared.html") }}, {{ make_url("truth estimated", "ckf_truth_estimated.html") }}
IVF: {{ make_url("seeded", "ivf_seeded.html") }}, {{ make_url("truth smeared", "ivf_truth_smeared.html") }}, {{ make_url("truth estimated", "ivf_truth_estimated.html") }}
{{ make_url("Ambiguity resolution", "ambi_seeded.html") }}
{{ make_url("Truth tracking", "truth_tracking.html") }}

### Vertexing {{ "" if all_exist(
    "vertexing_mu_scan.pdf", 
    "ivf_seeded_plots",
    "ivf_truth_smeared_plots",
    "ivf_truth_estimated_plots",
) else ":x: "}}

{% call detail_block("Vertexing vs. mu", "vertexing_mu_scan.pdf") %}
{{ make_image("vertexing_mu_scan.pdf", 350) }}
{% endcall %}

{% for mode in ["seeded", "truth_smeared", "truth_estimated"] %}

{% call detail_block("IVF "+mode, "ivf_"+mode+"_plots") %}
    
{% for url in [
    "covXX.pdf",
    "covYY.pdf",
    "diffx.pdf",
    "diffy.pdf",
    "diffz.pdf",
    "recoOverTrue.pdf",
] -%}
{{- make_image("ivf_"+mode+"_plots/"+url, "50%") -}}
{%- endfor %}

{% endcall %}

{% endfor %}


### CKF {{ "" if all_exist(
    "ckf_seeded_plots",
    "ckf_truth_smeared_plots",
    "ckf_truth_estimated_plots",
) else ":x: "}}

{% for mode in ["seeded", "truth_smeared", "truth_estimated"] %}

{% call detail_block("CKF "+mode, "ckf_"+mode+"_plots") %}
    
{% for url in [
    "trackeff_vs_eta.pdf",
    "trackeff_vs_pT.pdf",
    "nHoles_vs_eta.pdf",
    "nMeasurements_vs_eta.pdf",
] -%}
{{- make_image("ckf_"+mode+"_plots/"+url, "50%") -}}
{%- endfor %}

{% endcall %}

{% endfor %}

### Ambiguity resolution {{ "" if exists("ambi_seeded_plots") else ":x: "}}

{% call detail_block("seeded", "ambi_seeded_plots") %}
    
{% for url in [
    "trackeff_vs_eta.pdf",
    "trackeff_vs_pT.pdf",
    "nHoles_vs_eta.pdf",
    "nMeasurements_vs_eta.pdf",
] -%}
{{- make_image("ambi_seeded_plots/"+url, "50%") -}}
{%- endfor %}

{% endcall %}

### Truth tracking {{ "" if exists("truth_tracking_plots") else ":x: "}}

{% call detail_block("Truth tracking", "truth_tracking_plots") %}
    
{% for url in [
    "nHoles_vs_eta.pdf",
    "nMeasurements_vs_eta.pdf",
    "pull_d0.pdf",
    "pull_z0.pdf",
    "pull_theta.pdf",
    "pull_phi.pdf",
    "pull_qop.pdf",
    "pull_t.pdf",
] -%}
{{- make_image("truth_tracking_plots/"+url, "50%") -}}
{%- endfor %}

{% endcall %}
