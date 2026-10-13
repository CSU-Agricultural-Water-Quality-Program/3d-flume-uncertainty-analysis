from pathlib import Path
import csv, hashlib, json

ROOT = Path.cwd().resolve()
def checked(relative):
    path=(ROOT/relative).resolve()
    if not path.is_relative_to(ROOT): raise RuntimeError(f'Outside workspace: {path}')
    return path
def read(relative): return checked(relative).read_text(encoding='utf-8-sig')
def write(relative,text):
    path=checked(relative); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text,encoding='utf-8')
def move(old,new):
    src,dst=checked(old),checked(new)
    if not src.is_file() or dst.exists(): raise RuntimeError(f'Unsafe move: {src} -> {dst}')
    dst.parent.mkdir(parents=True,exist_ok=True)
    src.rename(dst)

# Preserve data and scientific outputs; hash them before any reorganization.
protected=[p for folder in ['data','output'] for p in checked(folder).rglob('*') if p.is_file()]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
write('tmp/meter_audit/reorganization_protected_before.json',json.dumps(manifest,indent=2))

original=read('code/bayesflume.r')
lines=original.splitlines()
locations=list(csv.DictReader(checked('tmp/meter_audit/reporting_function_locations.csv').open()))
header='''# Supporting reporting functions; sourced by code/base_flume.R.
# Defines summaries, diagnostics, rating-curve plots, and manuscript tables.
# Loading this file does not fit a model or export files.

'''
sections=[]
for row in locations:
    sections.append(f"# {row['name']}\n"+'\n'.join(lines[int(row['first'])-1:int(row['last'])]))
write('code/supporting_code/flume-reporting.R',header+'\n\n'.join(sections)+'\n')

# Remove definition regions, including their adjacent documentation, from entry.
regions=[(original.index('# Summaries use native draws'), original.index('# set working directory')),
         (original.index('######################### Graphing Outputs'),original.index('# Optional manufacturer coefficients')),
         (original.index('# Generative predictions output'),original.index('# grid and manufacturer equation')),
         (original.index('# Manuscript Table 4:'),original.index('table4 <- generate_manuscript_table4'))]
main=original
for start,end in sorted(regions,reverse=True): main=main[:start]+main[end:]
bootstrap='''# Resolve the project from this file, including when sourced in RStudio.
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.cli_file <- grep("^--file=", commandArgs(), value=TRUE)
.entry_file <- if (length(.source_files)) tail(unlist(.source_files),1L) else
  if (length(.cli_file)) sub("^--file=","",.cli_file[1]) else
    stop("Run or source code/base_flume.R as a file.")
flume_repo_dir <- dirname(dirname(normalizePath(.entry_file,winslash="/")))
source(file.path(flume_repo_dir,"code/supporting_code/project-paths.R"))
flume_repo_dir <- flume_project_root(flume_repo_dir)

# 1. Load fixed settings and supporting functions (no fitting here).
source(file.path(flume_repo_dir,"code/supporting_code/flume-metric-units.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-reporting.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-point-plots.R"))
source(file.path(flume_repo_dir,"code/supporting_code/flume-uncertainty-components.R"))
'''
start=main.index('flume_repo_dir <- if')
end=main.index('# Sampling requires an explicit opt-in',start)
main=main[:start]+bootstrap+main[end:]
main=main.replace('# bayesflume.r','# base_flume.R — main analysis entry point',1)
main=main.replace('# set working directory\nsetwd(file.path(flume_repo_dir, "data"))',
                  '# 2. Set input/output paths without changing the working directory.')
main=main.replace('output_dir <- "../output/bayes/oct2026_metric"',
                  'output_dir <- file.path(flume_repo_dir,"output/bayes/oct2026_metric")')
main=main.replace('input_csv <- "New Data Oct 2026/flume_data_oct2026_metric.csv"',
                  'input_csv <- file.path(flume_repo_dir,"data/New Data Oct 2026/flume_data_oct2026_metric.csv")')
main=main.replace('code/prepare-flume-metric-data.r','code/tools/prepare-flume-metric-data.R')
main=main.replace('# Define the Bayesian models','# 3. Define and fit the Bayesian model; archived alternatives stay commented.')
main=main.replace('source("../code/flume-point-plots.r")\n','')
main=main.replace('source("../code/flume-uncertainty-components.r")\n','')
main=main.replace('# Optional manufacturer coefficients by flume index (edit as needed)',
                  '# 4. Export rating-curve, point, and uncertainty-component figures.')
main=main.replace('# grid and manufacturer equation','# 5. Export stage-discharge predictions and manuscript Table 4.')
main=main.replace('table4 <- generate_manuscript_table4','# Parameter summaries use all saved draws.\ntable4 <- generate_manuscript_table4')
move('code/bayesflume.r','code/base_flume.R')
write('code/base_flume.R',main)

moves={
 'code/flume-metric-units.r':'code/supporting_code/flume-metric-units.R',
 'code/flume-point-plots.r':'code/supporting_code/flume-point-plots.R',
 'code/flume-uncertainty-components.r':'code/supporting_code/flume-uncertainty-components.R',
 'code/prepare-flume-metric-data.r':'code/tools/prepare-flume-metric-data.R',
 'code/validate-flume-metric.r':'code/tools/validate-flume-metric.R',
 'code/regenerate-flume-point-plots.r':'code/tools/regenerate-flume-figures.R',
 'code/preview-flume-uncertainty-components.r':'code/_unused_old_r/preview-flume-uncertainty-components.R',
}
for old,new in moves.items(): move(old,new)
python_files=list(checked('code').glob('*.py'))
for path in python_files:
    new='code/_unused_old_python/'+path.name
    move(str(path.relative_to(ROOT)),new)
    s=read(new).replace('Path(__file__).resolve().parents[1]','Path(__file__).resolve().parents[2]')
    if path.name=='main.py': s=s.replace('REPO = HERE.parent','REPO = HERE.parent.parent')
    s=s.replace('python code/simulate_flume_data_model3.py',
                'python code/_unused_old_python/simulate_flume_data_model3.py')
    s=s.replace('this file lives in <root>/code/.','this file lives in <root>/code/_unused_old_python/.')
    s=s.replace('bayesflume.r','base_flume.R (historical Model 3)')
    write(new,s)

replacements={**moves,'code/bayesflume.r':'code/base_flume.R'}
# Active documentation and code paths follow the new layout; audit narratives
# and source manifests remain historical provenance rather than being rewritten.
for name in ['README.md','output/README.md',*moves.values()]:
    s=read(name)
    for old,new in replacements.items(): s=s.replace(old,new)
    s=s.replace('bayesflume.r','base_flume.R')
    s=s.replace('code/main.py','code/_unused_old_python/main.py')
    for path in python_files:
        s=s.replace('code/'+path.name,'code/_unused_old_python/'+path.name)
    # Tools/archive entry files are one directory deeper than the former code root.
    if name.startswith('code/tools/') or name.startswith('code/_unused_old_r/'):
        s=s.replace('dirname(dirname(normalizePath(', 'dirname(dirname(dirname(normalizePath(')
        s=s.replace('winslash="/")))','winslash="/"))))')
        s=s.replace('winslash = "/")))','winslash = "/"))))')
    write(name,s)

# Regeneration imports reporting definitions directly rather than parsing entry.
s=read('code/tools/regenerate-flume-figures.R')
start=s.index('expressions <- parse(')
end=s.index('g <- plot_flume_posterior',start)
s=s[:start]+'source(file.path(repo_dir,"code/supporting_code/flume-reporting.R"))\n'+s[end:]
write('code/tools/regenerate-flume-figures.R',s)

# Validator imports extracted reporting functions; model/data still come from entry.
s=read('code/tools/validate-flume-metric.R')
start=s.index('# Only function definitions are loaded')
end=s.index('source("code/supporting_code/flume-point-plots.R")',start)
s=s[:start]+'# Import reporting functions without executing the main workflow.\nsource("code/supporting_code/flume-reporting.R")\n'+s[end:]
write('code/tools/validate-flume-metric.R',s)
print(f'Organized main entry, nine reporting functions, seven R helpers/tools, and {len(python_files)} obsolete Python scripts.')
