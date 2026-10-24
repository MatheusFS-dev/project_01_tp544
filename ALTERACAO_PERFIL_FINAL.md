# Ajuste do perfil final para tempo de execução

## Motivo

O perfil final anterior combinava uma malha local de 2 nm sobre toda a célula
periódica de 320 x 320 nm, uma extensão local de 400 nm em `z`, um domínio de
4 um em `z` e tempo máximo de 10 ps. Na GPU local, a estimativa exibida para o
primeiro caso (`A_air`, 0 graus) era de várias horas, inviável para a meta de
no máximo aproximadamente 15 minutos por simulação.

## Alterações no perfil `final`

| Item | Valor anterior | Valor atual | Justificativa |
| --- | ---: | ---: | --- |
| Passo da malha local | 2 nm | 3 nm | Mantém cerca de 13 células na largura de 40 nm e 10 na espessura de 30 nm do anel. |
| Região de malha fina em `x` e `y` | 320 x 320 nm | 220 x 220 nm | A malha fina cobre o anel de raio externo de 90 nm, com margem de 20 nm. |
| Região de malha fina em `z` | -200 a 200 nm | -30 a 60 nm | Mantém refinamento no metal e nas interfaces imediatas, sem aplicá-lo longe da estrutura. |
| Domínio FDTD em `z` | -2000 a 2000 nm | -800 a 800 nm | Conserva 280 nm entre o monitor de reflexão e a PML e 300 nm entre o monitor de transmissão e a PML. |
| Tempo máximo | 10 ps | 10 ps | O limite temporal é preservado para não interromper o decaimento residual prematuramente. |
| Auto-shutoff | `1e-6` | `1e-6` | A simulação só termina antecipadamente quando a energia residual tiver decaído ao critério definido. |
| PML | 16 camadas | 16 camadas | Valor mantido. |

O maior ganho vem de limitar a malha de 3 nm à vizinhança do anel. A região
refinada passa de toda a célula e 400 nm de altura para 220 x 220 x 90 nm,
reduzindo fortemente o número de células finas.

## GPU e fonte

- O perfil final executa com `run("FDTD","GPU")` nos espectros, na referência
  e nos campos próximos.
- BFAST não é usado.
- O modo `fast` permanece inalterado: malha local de 6 nm, domínio de -1200 a
  1200 nm, tempo máximo de 2,5 ps e execução `run;`.

## Como aplicar

1. Interromper qualquer execução criada com a configuração anterior usando
   **Quit & don't save**.
2. Remover os resultados finais antigos, preservando os arquivos `fast_`.
3. Executar `01_build_geometry;` para gerar novamente `A_air.fsp` e
   `B_SiO2.fsp` com a nova malha.
4. Executar `02_run_spectra;`.
5. Confirmar no log de cada caso que a energia residual atingiu o auto-shutoff;
   se uma execução chegar a 10 ps sem atingir o critério, esse resultado deve
   ser avaliado antes de ser aceito.

O tempo de execução depende da GPU e não pode ser garantido só pelo script,
mas esta configuração foi reduzida para tornar a meta de aproximadamente 15
minutos por simulação viável sem voltar ao perfil `fast`.

