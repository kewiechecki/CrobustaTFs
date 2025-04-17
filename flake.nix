{
  description = "Flake to build + develop CrobustaTFs R package";

  nixConfig = {
    bash-prompt = "\[CrobustaTFs$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    nixpkgs.url     = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs  = import nixpkgs { inherit system; config.allowUnfree = true; };
        rpkgs = pkgs.rPackages;
#		url = "https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FKH-ENS.blast.zip";
#		zipfile = "'index.html?file=data%2Fcirobu%2FKH-ENS.blast.zip'";
#		ensembl = "KH-ENS.blast";
#
#		downloadHook = ''
#		wget --no-check-certificate ${url}
#		unzip -o ${zipfile}
#		mv ${ensembl} data-raw
#		'';

        myPkg = rpkgs.buildRPackage rec {
          name    = "CrobustaTFs";
          version = "0.0.3";
          src     = ./.;

          # R‐side dependencies:
          propagatedBuildInputs = [
            rpkgs.devtools
            rpkgs.motifmatchr
            rpkgs.TFBSTools
            rpkgs.XML
          ];

          # Make sure R is on PATH at build time, and pkg-config can see libxml2 & gsl:
          nativeBuildInputs = [
            pkgs.R
            pkgs.pkg-config
          ];

          # C‑library dependencies for XML.so and DirichletMultinomial.so
          buildInputs = [
			pkgs.bzip2
			pkgs.curl
            pkgs.gsl
			pkgs.icu75
			pkgs.libpng
            pkgs.libxml2
          ];

		  preBuild = ''
		  make build
		  mv build $out
		  '';

          # re‑enable Nix’s R-wrapper so it injects R_LD_LIBRARY_PATH
          dontUseSetLibPath = false;

          meta = with pkgs.lib; {
            description = "…";
            license     = licenses.mit;
            maintainers = [ maintainers.kewiechecki ];
          };
        };
      in rec {
        # 1) allow `nix build` with no extra attr:
        defaultPackage = myPkg;

        # 2) drop you into a shell for interactive R work:
        devShells = {
          default = pkgs.mkShell {
            name = "crobusta-tfs-shell";
            buildInputs = [
              pkgs.git
              pkgs.R
              rpkgs.devtools
              rpkgs.motifmatchr
              rpkgs.TFBSTools
              rpkgs.XML
			  pkgs.bzip2
			  pkgs.curl
              pkgs.gsl
			  pkgs.icu75
			  pkgs.libpng
              pkgs.libxml2
              pkgs.pkg-config
            ];
            shellHook = ''
source ${pkgs.git}/share/bash-completion/completions/git-prompt.sh
# ensure that at _runtime_ R can find libxml2.so.2 and libgsl.so.28
export LD_LIBRARY_PATH="${pkgs.libpng.out}/lib:${pkgs.icu75.out}/lib:${pkgs.bzip2.out}/lib:${pkgs.curl.out}/lib:${pkgs.libxml2.out}/lib:${pkgs.gsl.out}/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="${pkgs.libpng.out}/lib:${pkgs.icu75.out}/lib:${pkgs.bzip2.out}/lib:${pkgs.curl.out}/lib/pkgconfig:${pkgs.libxml2.out}/lib/pkgconfig:${pkgs.gsl.out}/lib/pkgconfig:$PKG_CONFIG_PATH"
# ensure pkg-config still sees our .pc files
echo "→ LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "→ PKG_CONFIG_PATH=$PKG_CONFIG_PATH"
            '';
          };
        };
      }
    );
}

