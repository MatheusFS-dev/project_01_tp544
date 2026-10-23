# Simulação dos nanoanéis de Au

Passo a passo para simular a metasuperfície em ar e sobre SiO₂. Os detalhes físicos e numéricos estão em [SIMULATION_PARAMETERS.md](SIMULATION_PARAMETERS.md).

## Preparação

Você precisa do **Ansys Lumerical FDTD Solutions** para executar os arquivos `.lsf` e de **Python** para analisar e plotar os resultados. Na raiz do projeto, instale as dependências Python e crie as pastas de saída:

```powershell
python -m pip install -r requirements.txt
New-Item -ItemType Directory -Force results/geometry, results/spectra, results/resonances, results/nearfields
```

No Lumerical, defina `src/lumerical_scripts` como diretório de trabalho. Execute os scripts `.lsf` desse diretório, na ordem abaixo.

## Passo a passo

1. **Criar a geometria:** execute `01_build_geometry.lsf`. Confira os modelos `results/geometry/A_air.fsp` e `B_SiO2.fsp`.
2. **Verificar as referências sem Au:** em `run_reference.lsf`, deixe `case_index = 1` e execute. Depois altere para `case_index = 2` e execute novamente. Os TSVs `reference_*.tsv` ficam em `results/spectra`. Em 0°, espere aproximadamente `T = 1` no ar e `R = 0,0337` sobre SiO₂.
3. **Simular os espectros:** execute `02_run_spectra.lsf`. Ele roda os dois meios nos ângulos 0°, 15°, 30°, 45°, 60° e 75°. Cada caso concluído gera FSP, TSV, MAT e um arquivo `.done` em `results/spectra`. Uma nova execução pula os casos já concluídos.
4. **Plotar e escolher ressonâncias:** no terminal, execute os comandos abaixo. A primeira execução cria `results/resonances/resonances.csv`. Abra os gráficos T/R/A gerados em `results/resonances`, preencha esse CSV com as ressonâncias observadas e execute o mesmo comando novamente. Use `A_air` ou `B_SiO2` na coluna `case`, um dos seis ângulos na coluna `angle_deg`, um identificador de modo em `resonance_id` e um comprimento de onda inteiro entre 450 e 1600 nm em `wavelength_nm`.

   ```powershell
   cd src/lumerical_scripts
   python 03_review_resonances.py
   ```

5. **Simular campos próximos:** depois de revisar `results/resonances/nearfield_requests.tsv`, execute `04_run_nearfields.lsf` no Lumerical. Ele usa apenas as ressonâncias escolhidas em 0° e 45°.
6. **Plotar os campos:** no terminal, a partir de `src/lumerical_scripts`, execute:

   ```powershell
   python 04_plot_nearfields.py
   ```

Os mapas XY e XZ em PNG e PDF ficam em `results/nearfields`. A tabela de ressonâncias e o gráfico de evolução angular ficam em `results/resonances`.

## Testes locais

Na raiz do projeto, execute `python -m pytest -q`. Esses testes verificam o processamento com dados sintéticos. A validação física e a convergência exigem executar as simulações no Lumerical.
