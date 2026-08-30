"""Recompute and check the artifact hash(es) in this folder.
Usage:  python verify_artifact.py"""
import glob, hashlib, os, sys
here = os.path.dirname(os.path.abspath(__file__))
bad = ok = 0
for sf in glob.glob(os.path.join(here, '*.sha256')):
    want, name = open(sf).read().split()[:2]
    fp = os.path.join(here, name)
    got = hashlib.sha256(open(fp, 'rb').read()).hexdigest()
    good = (got == want)
    ok += good; bad += (not good)
    print(('MATCH   ' if good else 'MISMATCH'), name)
    if not good:
        print('  expected', want)
        print('  computed', got)
print('RESULT:', 'ALL MATCH' if bad == 0 else '%d MISMATCH' % bad)
sys.exit(0 if bad == 0 else 1)
