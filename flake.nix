{
  inputs = {
    nixpkgs.url = "github:cachix/devenv-nixpkgs/rolling";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
    devenv.inputs.nixpkgs.follows = "nixpkgs";
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, devenv, systems, ... } @ inputs:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      packages = forEachSystem (system: {
        devenv-up = self.devShells.${system}.default.config.procfileScript;
      });

      devShells = forEachSystem
        (system:
          let
            pkgs = nixpkgs.legacyPackages.${system};
          in
          {
            default = devenv.lib.mkShell {
              inherit inputs pkgs;
              modules = [
                {
                  # https://devenv.sh/reference/options/
                  packages = with pkgs; with rPackages; [
                    icu.dev
                    languageserver
                    glibcLocales
                    compositions
                    tidyverse
                    dotenv
                    data_table
                    mice
                    survival
                    rms
                    fastglm
                    cowplot
                    extrafont
                    stringr
                    styler
                    svglite
                  ];

                  languages.r = {
                    enable = true;
                  };
                  env.R_LIBS_USER="./.Rlib";

                  pre-commit.hooks.formatter = {
                    enable = true;
                    entry = ''
                    Rscript -e "styler::style_dir(path = '.',
                      recursive = TRUE,
                      filetype = c('R', 'Rmd'),
                      exclude_dirs = c('.devenv', '.direnv', '.git'),
                      transformers = styler::tidyverse_style())"'';
                  };
                }
              ];
            };
          });
    };
}
