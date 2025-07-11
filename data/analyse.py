import sys, os, urllib
try:
    import dune
except ImportError:
    dune = None

from importlib.metadata import version

def getVersion():
    femVersion = version("dune.fem")
    while True:
        dataPath = os.path.join(".",femVersion)
        if os.path.exists(dataPath):
            return femVersion
        ver = femVersion.split(".")
        if "dev" in ver[-1]:
            ver[-1] = 100 # try to find the largest matching release version
        lower = int(ver[-1])-1
        if lower < 0:
            # final try with only major.minor
            femVersion = femVersion.split(".")[0] + "." + femVersion.split(".")[1]
            dataPath = os.path.join(".",femVersion)
            if os.path.exists(dataPath):
                return femVersion
            # we give up
            raise RuntimeError(f"no suitable data for {version('dune.fem')} found")
        ver[-1] = str(lower)
        femVersion = ".".join(ver)

def main():
    femVersion = getVersion()
    dataPath = os.path.join(".",femVersion)
    print(f"obtaining data from {dataPath}")

    scripts = [m for m in os.listdir(dataPath) if m.endswith(".modules")]

    check = {}
    for script in scripts:
        filename = os.path.join(dataPath,script)
        with open(filename,"r") as f:
            modules = f.readlines()
        modules = set([c.strip()+".cc" for c in modules])
        scr = script.split(".modules")[0]
        for m in modules:
            if m in check:
                check[m] += [scr]
            else:
                check[m] = [scr]
    print("modules needed for 2 to 3 scripts:")
    for m in check:
        if len(check[m])==2 or len(check[m])==3:
            print(m,check[m])
    print("modules needed for 4 to 6 scripts:")
    for m in check:
        if len(check[m])>3 and len(check[m])<7:
            print(m,check[m])

if __name__ == "__main__":
    main()
