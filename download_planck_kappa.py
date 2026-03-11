#!/usr/bin/env python3
"""Download Planck 2018 CMB lensing convergence (kappa) map."""
import os, sys, subprocess

DATA_DIR = 'data'
os.makedirs(DATA_DIR, exist_ok=True)

# Try multiple sources for the Planck lensing map
urls = [
    # Planck 2018 (PR3) lensing - convergence map
    ('https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/lensing/COM_Lensing_4096_R3.00.tgz', 'planck_lensing_r3.tgz'),
    # Planck 2015 (PR2) lensing
    ('https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/lensing/COM_CompMap_Lensing_2048_R2.00.fits', 'planck_lensing_r2.fits'),
    # NASA LAMBDA mirror
    ('https://lambda.gsfc.nasa.gov/data/foregrounds/lensing/COM_Lensing_4096_R3.00.tgz', 'planck_lensing_lambda.tgz'),
    # ESA PLA direct
    ('http://pla.esac.esa.int/pla/aio/product-action?LENSING.FILE_ID=COM_Lensing_4096_R3.00.tgz', 'planck_lensing_esa.tgz'),
]

for url, fname in urls:
    fpath = os.path.join(DATA_DIR, fname)
    if os.path.exists(fpath) and os.path.getsize(fpath) > 10000:
        print(f"Already have {fpath} ({os.path.getsize(fpath)} bytes)")
        sys.exit(0)
    
    print(f"Trying: {url[:80]}...")
    result = subprocess.run(
        ['curl', '-sL', '--max-time', '120', '-o', fpath, url],
        capture_output=True, text=True
    )
    
    if os.path.exists(fpath):
        size = os.path.getsize(fpath)
        # Check if it's an error page
        if size < 10000:
            with open(fpath, 'rb') as f:
                header = f.read(100)
            if b'HTML' in header or b'html' in header or b'404' in header:
                print(f"  Got HTML error page ({size} bytes), removing")
                os.remove(fpath)
                continue
        print(f"  Downloaded {fpath} ({size} bytes)")
        
        # If it's a tarball, extract it
        if fname.endswith('.tgz'):
            extract_dir = os.path.join(DATA_DIR, 'planck_lensing')
            os.makedirs(extract_dir, exist_ok=True)
            print(f"  Extracting to {extract_dir}...")
            subprocess.run(['tar', 'xzf', fpath, '-C', extract_dir], check=True)
            print(f"  Extracted. Contents:")
            for root, dirs, files in os.walk(extract_dir):
                for f in files:
                    fp = os.path.join(root, f)
                    print(f"    {fp} ({os.path.getsize(fp)} bytes)")
        sys.exit(0)
    else:
        print(f"  Download failed")

print("All download attempts failed!")
print("Falling back to building proxy from SDSS foreground galaxy density...")
sys.exit(1)
