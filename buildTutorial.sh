# setup venv like this
# python3 -m venv venv
# pip install -U pip
# pip install -r ./reqCI.txt

cd /host

source $HOME/dune-env/bin/activate
python ~/DUNE/dune-common/bin/rmgenerated.py --all
cd doc

make clean
rm generatorCompiler.*
export DUNE_SAVE_BUILD=APPEND
make -j 2 -i >& tutorial.out
