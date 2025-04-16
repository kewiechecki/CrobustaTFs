{

  nixConfig = {
    bash-prompt = "\[CrobustaTFs$(__git_ps1 \" (%s)\")\]$ ";
  };

  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs, utils }: 
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };

      in {
        defaultPackage = pkgs.stdenv.mkDerivation {
          name = "CrobustaTFs";
          src = ./.;
          buildInputs = [];
        };

        devShell = with pkgs; mkShell {
          name = "crobusta-tfs-shell";
          buildInputs = [
            git      # for git prompt support
            R
            pkgs.rPackages.devtools
          ];
          shellHook = ''
source ${git}/share/bash-completion/completions/git-prompt.sh

make install
		  '';
		};
	});
}
