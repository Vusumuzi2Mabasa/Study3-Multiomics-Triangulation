"""Quarto Word report from the shared research manuscript."""
from pathlib import Path
import argparse,json,re,shutil,subprocess
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.style import WD_STYLE_TYPE

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--quarto',default=shutil.which('quarto'));a=ap.parse_args()
    if not a.quarto:ap.error('Quarto is required')
    work=ROOT/'work/word';work.mkdir(parents=True,exist_ok=True)
    reference=Document();section=reference.sections[0]
    section.page_width=Inches(8.5);section.page_height=Inches(11)
    section.top_margin=section.bottom_margin=Inches(.7);section.left_margin=section.right_margin=Inches(.75)
    for name in ['First Paragraph','Compact','Source Code']:
        if name not in reference.styles:reference.styles.add_style(name,WD_STYLE_TYPE.PARAGRAPH).base_style=reference.styles['Normal']
    for style in reference.styles:
        if style.type==1:
            style.font.name='Calibri';style.font.size=Pt(11);style.font.color.rgb=RGBColor(0,0,0)
            style.paragraph_format.space_after=Pt(6);style.paragraph_format.line_spacing=1.04
    for name,size in [('Title',23),('Subtitle',12),('Heading 1',16),('Heading 2',13),('Heading 3',11)]:
        reference.styles[name].font.size=Pt(size);reference.styles[name].paragraph_format.keep_with_next=True
    for e in list(reference.styles.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    reference.save(work/'reference.docx')
    readme=(ROOT/'README.md').read_text(encoding='utf-8');body=readme.split('<!-- manuscript -->',1)[1].lstrip()
    body=body.split('## Reproducibility',1)[0]+'\n## References'+readme.split('## References',1)[1]
    disclosure=next(line for line in readme.splitlines() if line.startswith('**Simulation study.**'))
    body=disclosure+'\n\n'+body;body=body.replace('](figures/','](../figures/')
    body=re.sub(r'!\[([^]]*)\]\(([^)]+)\)',r'![](\2){fig-alt="\1" width=6.4in}',body)
    repo='https://github.com/Vusumuzi2Mabasa/Study3-Multiomics-Triangulation/blob/main/'
    body=re.sub(r'(?<!!)\[([^]]+)\]\((?!https?://|#)([^)]+)\)',lambda m:f'[{m[1]}]({repo}{m[2]})',body)
    header='''---
title: "Triangulation of Genotype and RNA Signals in Simulated Breast Cancer"
author: "Vusumuzi Mabasa"
subtitle: "Complementary molecular discovery signatures and genetic associations"
format:
  docx:
    reference-doc: ../work/word/reference.docx
    toc: false
execute:
  enabled: false
---

'''
    qmd=ROOT/'report/manuscript.qmd';qmd.write_text(header+body,encoding='utf-8',newline='\n')
    subprocess.run([a.quarto,'render',qmd.name,'--to','docx','--output','Study3_Triangulation_Report.docx'],cwd=qmd.parent,check=True)
    output=qmd.parent/'Study3_Triangulation_Report.docx';doc=Document(output);doc.core_properties.author='Vusumuzi Mabasa'
    for style in doc.styles:
        if style.type==1:style.font.color.rgb=RGBColor(0,0,0)
    for p in doc.paragraphs:
        if p.style.name.startswith('Heading') or p._p.xpath('.//w:drawing') or p.text.startswith('Table '):p.paragraph_format.keep_with_next=True
        if p._p.xpath('.//w:drawing'):p.alignment=1;p.paragraph_format.space_before=Pt(10)
        if p.text.startswith('Figure '):p.paragraph_format.space_after=Pt(9)
    for shape in doc.inline_shapes:
        ratio=Inches(6.4)/shape.width;shape.width=Inches(6.4);shape.height=int(shape.height*ratio)
    for old in list(doc.tables):
        new=doc.add_table(rows=len(old.rows),cols=len(old.columns));new.style='Table Grid'
        for source,dest in zip(old.rows,new.rows):
            for a,b in zip(source.cells,dest.cells):b.text=a.text
        old._tbl.addprevious(new._tbl);old._tbl.getparent().remove(old._tbl)
    for table in doc.tables:
        table.autofit=True
        for i,row in enumerate(table.rows):
            row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
            if i==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
            for cell in row.cells:
                shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E8EFF4' if i==0 else ('F5F7F9' if i%2 else 'FFFFFF'));cell._tc.get_or_add_tcPr().append(shade)
                for p in cell.paragraphs:
                    p.paragraph_format.space_before=Pt(3);p.paragraph_format.space_after=Pt(3)
                    if len(table.rows)<=10:p.paragraph_format.keep_with_next=i<len(table.rows)-1
                    for r in p.runs:r.font.size=Pt(9);r.bold=i==0
    for style in doc.styles:
        for e in list(style.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    for section in doc.sections:
        p=section.footer.paragraphs[0];p.alignment=2;p.add_run('Vusumuzi Mabasa | Simulated data | ')
        field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');p._p.append(field)
        for r in p.runs:r.font.size=Pt(8)
    doc.save(output)
    structure=dict(author=doc.core_properties.author,figures=len(doc.inline_shapes),tables=len(doc.tables),editable_equations=len(doc.element.xpath('.//m:oMath')))
    (ROOT/'results/document_structure.json').write_text(json.dumps(structure,indent=2)+'\n',encoding='utf-8',newline='\n');print(structure)


if __name__=='__main__':main()
