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
    if dune:
        from dune.fem.utility import FemThreadPoolExecutor
        from dune.generator import builder

        modules = ["concepts"]
        workers = 4
        if len(sys.argv)>1:
            workers = int(sys.argv[1])
            if len(sys.argv)>2:
                modules = sys.argv[2:]

        builder.initialize()

        femVersion = getVersion()
        dataPath = os.path.join(".",femVersion)
        print(f"obtaining data from {dataPath}")

        onlyCompile = []
        if not "all" in modules:
            for m in modules:
                filename = os.path.join(dataPath,m+".modules")
                try:
                    with open(filename,"r") as f:
                        onlyCompile += f.readlines()
                except FileNotFoundError:
                    pass
            onlyCompile = set([c.strip()+".cc" for c in onlyCompile])
        else:
            onlyCompile = os.listdir(dataPath)

        with FemThreadPoolExecutor(max_workers=workers) as executor:
            for file in onlyCompile:
                filename = os.path.join(dataPath,os.fsdecode(file))
                if ( not filename.endswith(".cc") or
                     not os.path.exists(filename) ):
                    continue
                with open(filename,"r") as f:
                    executor.submit( builder.load, file.replace(".cc",""), f.read(), None )

        """
        # make sure things are safe on Mac
        # https://docs.python.org/3/library/urllib.request.html
        os.environ["no_proxy"] = "*"

        url = "https://gitlab.dune-project.org/dune-fem/ddfem/-/raw/improveDocBuild/doc/XXX?ref_type=heads"
        if module == "all":
            modules = ["boundary","examples","intro",
                       "transformers", "geometry"]
        else:
            modules = [module]
        for filebase in modules:
            a = urllib.request.urlretrieve(
                 url.replace("XXX",f"{filebase}.py"),
                 filename=f"{filebase}.py")
            b = urllib.request.urlretrieve(
                 url.replace("XXX",f"{filebase}_nb.ipynb"),
                 filename=f"{filebase}_nb.ipynb"
                 )
        """

if __name__ == "__main__":
    main()
