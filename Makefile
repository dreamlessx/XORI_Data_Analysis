# XORI cross-orientation analysis pipeline.
# All targets assume cwd = repo root and a populated .venv (Python 3.12).
# Run `make env` first if .venv is missing.

PY := .venv/bin/python
PIP := .venv/bin/pip

.PHONY: env figures stats supp manuscript all clean help

help:
	@echo "XORI pipeline targets:"
	@echo "  make env          create .venv and install requirements.txt"
	@echo "  make metrics      regenerate metric_data/all_roi from raw_data/tc_data"
	@echo "  make culls        regenerate metric_data/{cull_roi,r_cull_roi}"
	@echo "  make depth        regenerate depth_data/ scatter and ROI maps"
	@echo "  make covariates   regenerate data_*/ (baseline, halfwidth, lhi, osi, size, spatial)"
	@echo "  make supp         regenerate supplementary_analysis/outputs/"
	@echo "  make stats        print canonical stats to stdout (LaTeX-formatted)"
	@echo "  make figures      regenerate paper/figures/fig1, fig3..fig9 (fig2 is separate)"
	@echo "  make fig2         build paper/figures/fig2_model from Methods spec (when implemented)"
	@echo "  make manuscript   latexmk build of paper/manuscript.tex (requires basictex/mactex)"
	@echo "  make all          full pipeline: metrics → culls → depth → covariates → supp → stats → figures → manuscript"
	@echo "  make clean        remove latex build artifacts and __pycache__"

env:
	@command -v uv >/dev/null 2>&1 || { echo "install uv first: brew install uv"; exit 1; }
	uv venv --python 3.12 .venv
	uv pip install --python .venv/bin/python -r requirements.txt

metrics:
	$(PY) scripts/m_calc/all_metric.py

culls:
	$(PY) scripts/m_calc/cull_metric.py
	$(PY) scripts/m_calc/r_cull_metric.py

depth:
	$(PY) scripts/d_calc/depth.py

covariates:
	$(PY) scripts/bm_calc/baseline.py
	$(PY) scripts/bm_calc/halfwidth.py
	$(PY) scripts/bm_calc/lhi.py
	$(PY) scripts/bm_calc/osi.py
	$(PY) scripts/bm_calc/size.py
	$(PY) scripts/bm_calc/spatial.py

supp:
	$(PY) supplementary_analysis/scripts/run_all_analyses.py
	$(PY) supplementary_analysis/scripts/additional_analyses.py
	$(PY) supplementary_analysis/scripts/extended_analyses.py

stats:
	$(PY) paper/compute_stats.py

figures:
	$(PY) paper/make_figures.py

fig2:
	@if [ -f paper/make_fig2_model_schematic.py ]; then \
		$(PY) paper/make_fig2_model_schematic.py; \
	else \
		echo "paper/make_fig2_model_schematic.py not yet built. See README.md and docs/PIPELINE.md."; \
		exit 1; \
	fi

manuscript:
	@command -v latexmk >/dev/null 2>&1 || { echo "latexmk not found. Install: brew install --cask basictex"; exit 1; }
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=build manuscript.tex

all: metrics culls depth covariates supp stats figures manuscript

clean:
	rm -rf paper/build
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name .pytest_cache -prune -exec rm -rf {} +
