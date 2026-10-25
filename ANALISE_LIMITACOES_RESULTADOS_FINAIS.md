# Análise das limitações dos espectros finais

## Decisão de continuidade

Os espectros gerados pelo perfil `final` serão usados nas próximas etapas do
fluxo: revisão de ressonâncias (`03_review_resonances.py`) e mapas de campo
próximo (`04_run_nearfields.lsf`). Esta decisão é motivada pelo custo de uma
nova campanha de simulações: a execução atual demora tipicamente cerca de 30
minutos por caso e, em alguns casos registrados nos logs, demorou mais.

Consequentemente, os resultados oblíquos devem ser apresentados e
interpretados como **exploratórios**, com as limitações abaixo explicitadas no
relatório.

## Arquivos produzidos

Foram concluídos os 12 espectros exigidos, para os dois meios e os seis
ângulos:

- `A_air`: 0, 15, 30, 45, 60 e 75 graus;
- `B_SiO2`: 0, 15, 30, 45, 60 e 75 graus.

Cada arquivo TSV contém 1151 pontos entre 450 e 1600 nm, com transmissão,
reflexão e absorção calculada por `A = 1 - R - T`.

## Evidências observadas

### Encerramento temporal

Todos os casos atingiram o limite máximo de 10 ps antes de alcançar o critério
de auto-shutoff de `1e-6`. Exemplos extraídos dos logs:

| Caso | Auto-shutoff ao final | Tempo de parede |
| --- | ---: | ---: |
| `A_air`, 0 graus | `2,13e-4` | 23,9 min |
| `A_air`, 15 graus | `3,70e-3` | 36,0 min |
| `A_air`, 45 graus | `2,16e-2` | 36,0 min |
| `A_air`, 75 graus | `6,91e-2` | 64,9 min |
| `B_SiO2`, 0 graus | `3,93e-4` | 47,9 min |
| `B_SiO2`, 15 graus | `4,38e-3` | 72,5 min |
| `B_SiO2`, 45 graus | `1,64e-2` | 63,9 min |
| `B_SiO2`, 75 graus | `6,40e-2` | 38,2 min |

Portanto, ainda havia energia eletromagnética no domínio quando os monitores
foram avaliados. Esse residual pode introduzir oscilações e erros de
normalização nos espectros.

### Conservação de energia

Em uma simulação física estável, `R`, `T` e `A` devem permanecer próximos do
intervalo de 0 a 1, exceto por erro numérico muito pequeno. Os resultados em
30, 45, 60 e 75 graus mostram valores não físicos, incluindo valores negativos
e valores maiores que 1. O caso mais severo é `A_air`, 30 graus, que alcança
`A = -16,57`.

Em incidência normal, os espectros são qualitativamente mais coerentes, mas
ainda apresentam pequenas regiões não físicas: `A_air`, 0 graus, chega a
`A = -0,037`, e `B_SiO2`, 0 graus, chega a `A = -0,064`.

## Causas prováveis

### 1. Fonte de banda larga com fronteiras Bloch em incidência oblíqua

As fronteiras Bloch são necessárias para impor a fase transversal de uma onda
incidente oblíqua em uma célula periódica. Sem BFAST, a fonte
`Bloch/periodic` de banda larga usada na GPU não mantém exatamente o mesmo
ângulo físico para todos os comprimentos de onda. A fase Bloch corresponde a
um vetor de onda transversal, enquanto o módulo do vetor de onda muda com o
comprimento de onda.

Esse efeito cresce com o ângulo. Ele é a causa mais provável de os espectros
de 30 a 75 graus apresentarem os maiores erros de conservação de energia e
picos espúrios.

### 2. Campo residual ao atingir 10 ps

O auto-shutoff não foi alcançado em nenhum caso. O resíduo é pequeno nos casos
normais, mas cresce fortemente em incidência oblíqua. Ao finalizar antes de o
campo decair, as transformadas dos monitores ainda contêm transientes e
interferências, que alteram `R`, `T` e, consequentemente, `A`.

### 3. Sensibilidade maior em ângulos altos

Em ângulos altos, a componente transversal do vetor de onda aumenta e o
acoplamento com a célula periódica, a discretização espacial e a absorção nas
PMLs se tornam mais sensíveis. Assim, uma pequena imprecisão de fase, malha ou
decaimento temporal pode produzir um erro mais visível nos monitores de
potência.

### 4. Malha e geometria curva

A malha local de 3 nm resolve a largura do anel com aproximadamente 13 células
e a espessura com 10 células. Ela é uma escolha de compromisso para limitar o
tempo de execução, mas uma seção circular de ouro ainda é aproximada por uma
grade cartesiana. Esse efeito de discretização pode deslocar ou alargar
ressonâncias, principalmente quando combinado com incidência oblíqua.

## Interpretação recomendada

- Usar os resultados de 0 graus como os mais confiáveis do conjunto atual.
- Tratar picos isolados em 15 graus ou mais como candidatos, não como prova
  definitiva de uma ressonância física.
- Não interpretar absorção negativa ou maior que 1 como absorção real; esses
  valores são indicadores numéricos.
- Nos gráficos do relatório, manter a definição `A = 1 - R - T`, mas registrar
  esta limitação no texto e evitar atribuir significado físico a regiões fora
  de `[0, 1]`.
- Identificar as ressonâncias para os mapas de campo pela feição espectral mais
  estável e visualmente observável, priorizando 0 e 45 graus conforme exigido
  pela tarefa.

## Próximas etapas aprovadas

1. Executar `python src/lumerical_scripts/03_review_resonances.py`.
2. Preencher `results/resonances/resonances.csv` com as ressonâncias escolhidas.
3. Executar `04_run_nearfields.lsf` para as escolhas em 0 e 45 graus.
4. Executar `python src/lumerical_scripts/04_plot_nearfields.py`.
5. Ao redigir o relatório, citar esta limitação numérica e distinguir picos
   estáveis de artefatos evidentes.

