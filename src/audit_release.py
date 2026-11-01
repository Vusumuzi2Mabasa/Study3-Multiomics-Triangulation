"""Read-only public release audit, suitable for GitHub CI."""
from pathlib import Path
import hashlib,json,re
from docx import Document

ROOT=Path(__file__).resolve().parents[1]


def main():
    count=0
    for line in (ROOT/'checksums.sha256').read_text(encoding='utf-8-sig').splitlines():
        digest,name=line.split('  ',1);assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name;count+=1
    audit=json.loads((ROOT/'results/audit.json').read_text());assert audit['status']=='passed' and audit['training_test_pairs']==200 and audit['distinct_cohort_seeds']==2400
    assert audit['MOFA_final_converged']==200
    assert json.loads((ROOT/'results/R_validation.json').read_text())['status']=='passed'
    qa=json.loads((ROOT/'results/report_qa.json').read_text());assert qa['status']=='passed'
    doc=Document(ROOT/'report/Study3_Triangulation_Report.docx');assert doc.core_properties.author=='Vusumuzi Mabasa' and len(doc.inline_shapes)==8 and len(doc.tables)==5
    assert len(doc.element.xpath('.//m:oMath'))>=10
    text=(ROOT/'README.md').read_text(encoding='utf-8');assert 'Simulation study.' in text and 'No patient data' in text
    for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',text):
        if not link.startswith(('https://','http://','#')):assert (ROOT/link).exists(),link
    assert len(list((ROOT/'results/final').glob('*.json')))==200
    assert len(list((ROOT/'results/calibration').glob('*.json')))==2000
    print(f'Release passed: {count} file hashes, 200 pairs, 2000 calibration cohorts, eight figures and five tables')


if __name__=='__main__':main()
