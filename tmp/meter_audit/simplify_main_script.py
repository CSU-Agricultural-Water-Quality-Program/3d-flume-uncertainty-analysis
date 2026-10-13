from pathlib import Path
import hashlib
import re

root = Path(r'C:/Users/ansle/OneDrive/Documents/GitHub/3d-flume-uncertainty-analysis').resolve()
old = root / 'code/base_flume.R'
archive = root / 'code/_unused_old_r/base_flume.R'
new = root / 'code/bayes_flume.R'
for path in (old, archive, new):
    assert path.resolve().is_relative_to(root)
assert old.is_file() and not archive.exists() and not new.exists()
old_bytes = old.read_bytes()
text = old_bytes.decode('utf-8-sig').replace('\r\n', '\n')

# Keep the active model and reporting; remove the three commented alternatives.
begin = text.index('# 3. Define and fit the Bayesian model;')
end = text.index('# Conditional stage-error model:', begin)
text = text[:begin] + '# 3. Fit the Bayesian calibration model.\n\n' + text[end:]
begin = text.index('# Model 4 is disabled and retained for reference.')
end = text.index('# 4. Export rating-curve, point, and uncertainty-component figures.', begin)
text = text[:begin] + text[end:]
text = text.replace('# base_flume.R — main analysis entry point',
                    '# bayes_flume.R — main calibration and reporting script')
text = text.replace('# Sampling requires an explicit opt-in; syntax/data checks do not sample.\n'
                    'RUN_CALIBRATION <- identical(Sys.getenv("FLUME_RUN_CALIBRATION"), "true")\n', '')
text = text.replace('if (!RUN_CALIBRATION) {\n'
                    '  stop("Metric inputs and model are prepared. Sampling is disabled; explicitly set FLUME_RUN_CALIBRATION=true to fit.")\n'
                    '}\n', '')
text = text.replace('# Current calibration products are separate from older Bayes and Python outputs.\n'
                    '# output_dir <- "../output/bayes/sep2026"\n'
                    '# output_prefix <- "sep2026"\n'
                    '# Keep the previous October fit intact while other uncertainty terms are reviewed.\n',
                    '# Save the metric calibration separately from historical fits.\n')
text = text.replace('# the discharge likelihood through H_true. Other priors remain under review.',
                    '# the discharge likelihood through H_true.')
text = text.replace('    # Previous September 2026 proportional-error prior, retained for reference:\n'
                    '    # cv_Q ~ dlnorm( log(2e-4) , 0.5 ),        # median 0.02% of reading\n'
                    '    # October revision: median 0.2% = 0.002; prior width remains 0.5.\n',
                    '    # Proportional-SD prior: median 0.2% = 0.002; log-scale SD 0.5.\n')

# Use unnumbered analysis objects; retain existing result filenames for helpers.
text = re.sub(r'\bflume_model_3\b(?!["_])', 'flume_model', text)
text = re.sub(r'\bp3\b', 'model_summary', text)
text = re.sub(r'\bpost3\b', 'post', text)

# Compact path setup still supports Source, Rscript, and interactive execution.
begin = text.index('# Resolve from this file when sourced, or from getwd() when run line by line.')
end = text.index('# 1. Load fixed settings and supporting functions (no fitting here).', begin)
bootstrap = '''# Locate the project from this file, or from the working directory for Run.
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.cli_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
.entry_file <- if (length(.source_files)) tail(unlist(.source_files), 1L) else
  if (length(.cli_file)) .cli_file[1] else NULL
flume_repo_dir <- if (is.null(.entry_file)) getwd() else
  dirname(normalizePath(.entry_file, winslash="/"))
while (!file.exists(file.path(flume_repo_dir, "code/supporting_code/project-paths.R"))) {
  .parent <- dirname(flume_repo_dir)
  if (identical(.parent, flume_repo_dir))
    stop("Set the working directory to the flume repository or source code/bayes_flume.R.")
  flume_repo_dir <- .parent
}
source(file.path(flume_repo_dir, "code/supporting_code/project-paths.R"))
flume_repo_dir <- flume_project_root(flume_repo_dir)

'''
text = text[:begin] + bootstrap + text[end:]
text = re.sub(r'\n{3,}', '\n\n', text)

# Verify these edits do not alter the active likelihood, priors, or sampling call.
def active_model(content, name):
    start = content.index(name + ' <- ulam(')
    end = content.index('\n# Preserve complete draws', start)
    block = content[start:end]
    return '\n'.join(line.strip() for line in block.splitlines()
                     if line.strip() and not line.lstrip().startswith('#'))
assert active_model(text, 'flume_model') == active_model(old_bytes.decode('utf-8-sig'), 'flume_model_3').replace('flume_model_3 <-', 'flume_model <-')
assert 'FLUME_RUN_CALIBRATION' not in text and 'RUN_CALIBRATION' not in text
assert not re.search(r'flume_model_[124]|# Model [124]', text)
assert all(name in text for name in ('save_flume_point_plots(', 'save_flume_uncertainty_components(', 'generate_manuscript_table4('))
old.rename(archive)
assert archive.read_bytes() == old_bytes
new.write_text(text, encoding='utf-8', newline='\n')
print(f'Archived {len(old_bytes.decode("utf-8-sig").splitlines())} lines; new main has {len(text.splitlines())} lines.')
print('Active model expression unchanged apart from its object name; archive preserved byte for byte.')
print('Archive SHA256: ' + hashlib.sha256(old_bytes).hexdigest())
