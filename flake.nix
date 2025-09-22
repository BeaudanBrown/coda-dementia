{
  outputs =
    {
    self,
    nixpkgs,
    nixpkgsUnstable,
    flake-utils,
    devshell,
    pre-commit-hooks,
    ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        pkgsUnstable = nixpkgsUnstable.legacyPackages.${system};
      in
        {
        checks = {
          pre-commit-check = pre-commit-hooks.lib.${system}.run {
            src = ./.;
            hooks = {
              air-fmt = {
                enable = true;
                entry = "${pkgs.air-formatter}/bin/air format";
                files = ".*\.[rR]$";
              };
            };
          };
        };
        devShells.default =
          let
            pkgs = nixpkgs.legacyPackages.${system} // {
              overlays = [ devshell.overlays.default ];
            };
          in
            pkgs.mkShell {
              inherit (self.checks.${system}.pre-commit-check) shellHook;
              env.R_LIBS_USER = "./.Rlib";
              buildInputs = [
                pkgs.bashInteractive
                self.checks.${system}.pre-commit-check.enabledPackages
              ];
              packages = with pkgsUnstable;
                [
                  R
                  quarto
                  just
                  air-formatter
                ]
                ++ (with pkgsUnstable.rPackages; [
                  data_table
                  languageserver
                  dotenv
                  compositions
                  tidyverse
                  mice
                  survival
                  rms
                  fastglm
                  cowplot
                  extrafont
                  stringr
                  styler
                  svglite
                  parallel
                  boot
                  rlang
                  targets
                  tarchetypes
                  usethis
                  qs2
                  emmeans
                  lme4
                  lubridate
                  tidyverse
                  mvtnorm
                  renv
                  extrafont
                  crew
                  parallelly
                  future
                  furrr
                  broom
                  RhpcBLASctl
                  autometric
                  ks
                  patchwork
                  ggrepel
                ]);
            };
      }
    );

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    nixpkgsUnstable.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    pre-commit-hooks.url = "github:cachix/git-hooks.nix";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
    devshell.url = "github:numtide/devshell";
  };
}
