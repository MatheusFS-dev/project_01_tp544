# Investigação e execução HFSS corrected

Data e hora: 2026-11-09T07:38:56-03:00

## Estado

**Campanha incompleta e pausada para decisão sobre a tabela de modos.** Foram exportados seis dos doze espectros. Quatro não têm a lacuna de modos identificada. Dois são diagnósticos, com a faixa limitada assinalada nos gráficos. Os seis restantes não foram executados.

| Caso | Ângulos | Situação |
| --- | --- | --- |
| Suspensa no ar | 0°, 15° | CSV, PNG, PDF e projetos resolvidos. Tabelas contêm os modos propagantes. |
| Suspensa no ar | 30°, 45° | Apenas diagnóstico. Falta TM(-1,0) abaixo de 480 e 546,274 nm, respectivamente. |
| Suspensa no ar | 60°, 75° | Pendente da decisão sobre a tabela. |
| Sobre SiO₂ | 0°, 15° | CSV, PNG, PDF e projetos resolvidos. Tabelas contêm os modos propagantes. |
| Sobre SiO₂ | 30°, 45°, 60°, 75° | Pendente da decisão sobre a tabela. |

A solicitação de autorização é específica: adaptar a tabela de modos às ordens propagantes dos ângulos pedidos, mantendo número de modos, material, malha, tolerâncias, setup e varredura do colega. Nenhuma adaptação dessa tabela foi realizada sem resposta.

A cobertura da tabela e a convergência adaptativa não são uma validação completa de T/R/A. Não foram realizados estudos independentes de refinamento, controles de Fresnel ou comparação com dissipação volumétrica. A identidade R + T + A = 1 é imposta pela definição de A.


## O que foi preservado

A implementação nova está em `src/hfss_corrected`. Os projetos, CSV e gráficos ficam em `results/hfss_corrected`. O fluxo antigo e seus resultados foram preservados.

O ponto de partida é `Downloads/Treinar/HFSS_Base_Dados.py`, geometria `G0246`. O anel já tem os raios 90/50 nm, altura 30 nm e período 320 nm pedidos no enunciado. As funções de material, portas, setup, varredura e extração espectral foram copiadas sem alterações, conferidas pela árvore sintática Python.

Foram mantidos Au Johnson–Christy, adaptação em 666,2 THz, máximo de 15 passagens, ΔS de 0,01, varredura interpolada com 401 pontos e tolerância de 0,5%, e 4/10 modos no ar/SiO₂. A seleção de incidência é TM00, em vez da média das duas polarizações do colega. SiO₂ tem índice 1,45 e termina diretamente na porta inferior, sem interface inferior com ar.

**Diferença de material:** o professor sugere Palik consistente entre os trabalhos. Esta reprodução conserva Au_JC do colega em atendimento à instrução de não modificar suas configurações. Ela não é uma comparação controlada com o Palik antigo.

## Primeiro resultado medido

Em 798,746 nm, o novo caso suspenso a 0° tem T = 0,053520, R = 0,672657 e A = 0,273823. O HFSS antigo tinha T = 0,786581, R = 0,063367 e A = 0,150052 em 798,780 nm. Os comprimentos de onda são próximos, mas diferentes.

O erro RMS entre os novos T/R/A e G0246 é 0,002776 após alinhar os comprimentos de onda por interpolação. O erro RMS em relação ao HFSS antigo é 0,242348. O novo espectro não reproduz a curva antiga. A banda larga próxima de 800 nm reapareceu. Os gráficos mostram os dados nativos sem suavização nem limitação artificial dos valores. A faixa nativa é 450,0037–1599,7463 nm, devido aos limites arredondados de 187,4–666,2 THz mantidos do colega.

O AEDT registrou convergência em duas passagens, com ΔS = 0,003469 e 4624 elementos. Isso demonstra o atendimento ao critério do colega. Não substitui um estudo independente de convergência espectral.

## Diferenças e problemas encontrados

1. O construtor antigo não remove o anel do volume de ar. O construtor do colega faz essa subtração. Na execução nova, os volumes medidos no AEDT concordam com os volumes analíticos, incluindo o ar com o anel removido. O ouro é resolvido internamente, verificado diretamente na propriedade do editor e no projeto salvo.
2. Foi construída uma cópia isolada usando o construtor antigo sem alterações e a tabela Palik existente. O AEDT retornou efetivamente `Solve Inside = false` para o ouro e volume de ar de 81.920.000 nm³, igual ao bloco completo sem retirar o anel de 527.787,566 nm³. Essa reprodução confirma os dois defeitos no construtor atual. Ela não é uma inspeção do projeto histórico resolvido, que não está no repositório.
3. A classificação de modos por um teste de Gamma existe no fluxo antigo e não no cálculo do colega. É uma diferença de extração que ainda não foi isolada como causa dos dentes.
4. O modelo do colega é uma referência a 0° no ar. A tabela automática de quatro modos contém TE00, TM00, TE(-1,0) e TE(0,-1). Ela omite TM(-1,0), que pode propagar no ar a 30° ou mais na parte curta da faixa. Preservar essa tabela limita a validade dos resultados oblíquos. Foi solicitada autorização específica para adaptar somente a tabela, conservando número de modos e demais parâmetros.

A reconstrução recuperou a resposta do colega. Como material, setup, domínio e extração também diferem do modelo antigo, ela não estabelece que apenas a sobreposição causou toda a discrepância. O mecanismo dos dentes antigos permanece sem isolamento causal. A ausência de TM(-1,0) está comprovada nas novas tabelas oblíquas preservadas do colega. A coincidência entre os dentes novos e essa faixa é evidência diagnóstica, não uma comparação antes/depois de corrigir a tabela.

A execução oblíqua foi pausada após esse diagnóstico. O driver verifica a tabela salva e recusa resolver novos casos com modos propagantes ausentes, sem alterar automaticamente a configuração solicitada.

## Arquivos de evidência

- `results/hfss_corrected/audits/original_builder/audit.json`: reprodução isolada do construtor antigo, incluindo Solve Inside falso e volumes sobrepostos.
- `results/hfss_corrected/audits/mode_coverage.json`: ordens propagantes ausentes nas tabelas salvas.
- `results/hfss_corrected/audits/A_air_0deg.json`: volumes, configuração, validação e comparação.
- `results/hfss_corrected/audits/A_air_0deg_convergence.csv`: convergência registrada pelo solver.
- `results/hfss_corrected/exports/G0246_comparison.csv`: curvas novas, do colega e antigas alinhadas.
- `results/hfss_corrected/projects/A_air_0deg.aedt`: projeto resolvido e sua base de resultados.

## Referências

- Enunciado: `docs/problem.pdf`.
- Plano: `docs/Plano_Investigacao_HFSS_Atualizado.md`.
- Referência recebida: `~/Downloads/Treinar/HFSS_Base_Dados.py`, SHA-256 `6134b6b5ef6162eac16af5809fdc452e1104f2fc44d89d76647d194523e870a6`.
- [Documentação Ansys sobre modos e inclusão dos modos propagantes](https://ansyshelp.ansys.com/public/views/secured/electronics/v252/en/subsystems/hfss/Content/HFSS/Modes.htm).

## Consulta aos resultados

- `results/hfss_corrected/plots/spectra_0_15deg.pdf` e `.png`: comparação dos quatro espectros concluídos sem lacuna de modos propagantes.
- `results/hfss_corrected/plots/comparison_A_air_0deg.pdf` e `.png`: novo HFSS, G0246 e HFSS antigo.
- `results/hfss_corrected/plots/A_air_30deg.*` e `A_air_45deg.*`: diagnósticos com faixa limitada sombreada.
- `results/hfss_corrected/exports/campaign_status.csv`: estado de todos os doze casos.
- Os CSV individuais preservam os 401 pontos nativos, sem suavização.

## Verificação executada

Os 48 testes existentes passaram. Quatro verificações temporárias adicionais passaram durante o desenvolvimento: construção dos dois meios sem sobreposição e aceitação/rejeição de tabelas de modos em casos de controle. Os testes e artefatos temporários foram removidos após a verificação. Os scripts novos compilam. As funções copiadas de material, portas, setup, varredura e extração mantêm a mesma árvore sintática das funções do colega.

Para retomar o código: `uv sync`, ativar `.venv/bin/activate` e executar `python3 src/hfss_corrected/run_campaign.py`. Enquanto a tabela não for adaptada, o driver recusa os casos com modos propagantes ausentes. Os casos já exportados sem essa lacuna são reutilizados.
