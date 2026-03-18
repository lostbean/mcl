{
  description = "mcl - Haskell library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, treefmt-nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        treefmtWrapper = treefmt-nix.lib.evalModule pkgs ./treefmt.nix;

        haskellStack = pkgs.haskellPackages.ghcWithPackages (ps: with ps; [
          cabal-install
          haskell-language-server
          fourmolu
          cabal-fmt
        ]);

      in
      {
        formatter = treefmtWrapper.config.build.wrapper;

        checks = {
          formatting = treefmtWrapper.config.build.check self;
        };

        devShells.default = pkgs.mkShell {
          buildInputs = [
            haskellStack
            pkgs.pkg-config
            pkgs.zlib
            treefmtWrapper.config.build.wrapper
            pkgs.nixpkgs-fmt
            pkgs.lefthook
          ];

          shellHook = ''
            echo "mcl Dev Environment Loaded"
            echo "GHC version: $(ghc --version)"
            echo "Lefthook version: $(lefthook version)"
            # Ensure lefthook is installed in the git repo
            if [ -d .git ] && [ ! -f .git/hooks/pre-commit ]; then
               lefthook install
            fi
          '';
        };
      }
    );
}
