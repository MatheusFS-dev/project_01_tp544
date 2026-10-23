# Parâmetros e execução da simulação

## Sistema físico

A célula periódica 3D tem 320 × 320 nm. O anel de Au tem raio externo de 90 nm, raio interno de 50 nm e altura de 30 nm, de z = 0 a 30 nm. Usa-se **Au (Gold) - Palik** nos dois casos, sem camada de adesão. O caso A tem ar dos dois lados. O caso B tem ar acima e SiO₂ de índice 1,45 abaixo de z = 0. O bloco de SiO₂ chega a z = −3000 nm, atravessando a PML inferior sem criar outra interface óptica.

A onda plana incide do ar em direção a −z, no plano xz, com polarização TM/p (`polarization angle = 0`). São simulados θ = 0, 15, 30, 45, 60 e 75°, medidos da normal. Em 0°, E é paralelo a x. Há periodicidade em x e y. A incidência normal usa fonte `Bloch/periodic`; os ângulos oblíquos usam BFAST para manter o ângulo ao longo da banda.

## Domínio e custo

O FDTD ocupa z = −2000 a +2000 nm, com PML de 16 camadas em z. A fonte fica em +300 nm, o monitor de reflexão em +520 nm e o de transmissão em −500 nm. A malha automática tem `accuracy = 2`; uma caixa local de z = −200 a +200 nm e de uma célula em x e y impõe passo de 4 nm nos três eixos. Essa caixa pode consumir muita memória. A duração máxima é 2500 fs e o limiar de desligamento automático é 10⁻⁶. Deve-se confirmar que o desligamento automático ocorreu e que BFAST permaneceu estável, sobretudo em 60° e 75°.

Os 2000 nm de distância da origem a cada limite equivalem a **1,25 × 1600 nm**, e não a 1,5 vezes o maior comprimento de onda. Essa escolha é um compromisso numérico, sujeito à comparação com ±2400 nm. A malha de 4 nm também precisa ser comparada com 3 nm. A orientação de [malha não uniforme da Ansys](https://optics.ansys.com/hc/en-us/articles/360034382634-Understanding-the-non-uniform-mesh-in-FDTD) embasa o uso da caixa de malha, mas não demonstra convergência deste modelo.

## Dados e normalização

Os espectros usam 1151 pontos de 450 a 1600 nm, com passo de 1 nm. Os monitores espectrais armazenam apenas potência integrada. Para incidência descendente, `T = -transmission(T_monitor)`, `R = transmission(R_monitor)` e `A = 1 − R − T`. Valores fora dos limites físicos permanecem nos arquivos para expor erros numéricos. Cada TSV tem cabeçalho `wavelength_nm`, `T`, `R`, `A`, seguido de 1151 linhas. O MAT legado contém as mesmas grandezas. O marcador `.done` indica que FSP, TSV e MAT foram escritos.

As ressonâncias são selecionadas **manualmente** nos espectros. O arquivo `results/resonances/resonances.csv` tem as colunas `case,angle_deg,resonance_id,wavelength_nm`; `case` vale `A_air` ou `B_SiO2`. A identificação de modos e sua correspondência entre os meios são responsabilidade da análise física. O script Python valida a grade, produz o resumo angular e solicita campos apenas para 0° e 45°.

Cada simulação de campo usa uma faixa de fonte de ±1 nm ao redor da ressonância e um ponto de monitor. Os mapas são o plano XY em z = +35 nm e o corte XZ em y = 0, de z = −200 a +200 nm. O MAT contém coordenadas em nm, `E2_xy`, `E2_xz`, e comprimentos de onda solicitado e amostrado. Para fonte de amplitude 1 V/m e normalização CW padrão, `E2` representa `|E|²/|E₀|²`. Confirme essa normalização com referência homogênea antes de quantificar máximos, especialmente em BFAST. Os gráficos compartilham escala de cores entre A e B para o mesmo identificador e ângulo quando ambos existem.

## Sequência

Execute a partir de `src/lumerical_scripts` no FDTD Solutions. Antes, crie as pastas `results/geometry`, `results/spectra`, `results/resonances` e `results/nearfields` na raiz do repositório. Essas pastas são ignoradas pelo Git.

1. `01_build_geometry.lsf` cria os dois modelos-base, sem fontes ou monitores.
2. `run_reference.lsf`, primeiro com `case_index = 1` e depois 2, verifica as referências sem Au em 0°.
3. `02_run_spectra.lsf` executa as 12 combinações e grava cada resultado imediatamente. Pode ser retomado, pois ignora casos com `.done`.
4. `python 03_review_resonances.py` cria `resonances.csv` vazio na primeira execução e plota espectros concluídos. Preencha o CSV manualmente e execute novamente. O arquivo manual não é sobrescrito. Revise `resonance_summary.csv`, `resonance_vs_angle` e `nearfield_requests.tsv`.
5. `04_run_nearfields.lsf` executa as solicitações numéricas. Depois, `python 04_plot_nearfields.py` cria PNG e PDF. Não altere `resonances.csv` enquanto houver solicitações em execução.

Execute `python -m pytest -q` para os testes locais. Sem o executável Lumerical, a geometria, os sinais e os resultados numéricos ainda exigem verificação no simulador.

## Verificações antes do relatório

Na referência sem Au, espere aproximadamente T = 1, R = 0 e A = 0 no ar. Na interface ar/SiO₂ em 0°, espere R ≈ 0,0337 e T ≈ 0,9663. Confira 1151 linhas, extremidades de 450 e 1600 nm e `T+R+A=1` dentro da precisão numérica. Registre diferenças espectrais e deslocamentos das ressonâncias ao comparar limites ±2000/±2400 nm e malhas 4/3 nm. Verifique também a sensibilidade da PML e a duração efetiva da simulação nos ângulos altos. Só mapas das seleções em 0° e 45° devem ser produzidos.
