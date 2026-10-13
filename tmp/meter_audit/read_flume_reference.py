from pathlib import Path
from pypdf import PdfReader

p = Path("D:/OneDrive - Colostate/Cabot,Perry's files - 3D Printing Flumes/REFERENCES/large_60-degree-v_trapezoidal_flume_discharge_table.pdf")
r = PdfReader(p)
print('PAGES',len(r.pages))
for i, page in enumerate(r.pages):
    print(f'\nPAGE {i+1}\n{page.extract_text()}')
