import sys, os, urllib
try:
    import dune
except ImportError:
    dune = None

from importlib.metadata import version

def getVersion():
    femVersion = version("dune.fem")
    tried = [femVersion]
    dataPath = os.path.join(".",femVersion)
    if os.path.exists(dataPath):
        return femVersion
    # femVersion might be 2.12.1.3 but jit modules are only stored for 2.12.0
    # Currently we just try 2.12.0 and 2.12 and then give up
    femVersion = femVersion.split(".")[:-1]
    femVersion = '.'.join(femVersion)
    tried += [femVersion]
    dataPath = os.path.join(".",femVersion)
    if os.path.exists(dataPath):
        return femVersion
    femVersion = femVersion.split(".")[:-1]
    femVersion = '.'.join(femVersion)
    tried += [femVersion]
    dataPath = os.path.join(".",femVersion)
    if os.path.exists(dataPath):
        return femVersion
    raise RuntimeError(f"no suitable data for {version('dune.fem')} found. Tried {tried}.")
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

        try:
            femVersion = getVersion()
        except Exception as e:
            print("no version found so no precompilation available")
            print(e)
            sys.exit(0)

        dataPath = os.path.join(".",femVersion)
        print(f"obtaining data from {dataPath}")

        onlyCompile = []
        if not "all" in modules:
            for m in modules:
                filename = os.path.join(dataPath,m+".jitmodules")
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
