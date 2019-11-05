```shell
cmsrel CMSSW_11_0_0_pre11_ROOT618
cd CMSSW_11_0_0_pre11_ROOT618
cmsenv
virtualenv PYTHON
source PYTHON/bin/activate
pip install pyspark
git clone https://github.com/JavierCVilla/PyRDF.git
cd PyRDF/
pip install .
cd ..
rm -rf PyRDF/
cd ..
time python ./TtSkim.py
```
