% =========================================================================
% PBL1 - Distribuição de Energia Elétrica - IFES Guarapari - 2026/1
% Grupo 1
%
% OBJETIVO 4: Fluxo de Potência - Método Backward-Forward Sweep (BFS)
% Modelos de Carga: Tipo Z (impedância constante),
%                   Tipo I (corrente constante),
%                   Tipo P (potência constante)
%
% Uso: execute este script diretamente no MATLAB.
%      Os resultados são impressos no Command Window e gerados
%      gráficos comparativos ao final.
% =========================================================================

clc; clear; close all;

% =========================================================================
% 1. PARÂMETROS DO SISTEMA
% =========================================================================
VLL_kV  = 11.4;                        % Tensão de linha [kV]
Vbase   = (VLL_kV / sqrt(3)) * 1e3;   % Tensão de fase base [V]

% Fasores da tensão nominal na subestação (sequência positiva, ABC)
% Va = Vbase /0°,  Vb = Vbase /-120°,  Vc = Vbase /+120°
alpha = exp(1j * 2*pi/3);
V_sub = Vbase .* [1+0j ; alpha^2 ; alpha];   % vetor coluna [3x1]

% =========================================================================
% 2. MATRIZ DE IMPEDÂNCIA DE FASE [Zabc] - por milha
%    Calculada via Equações de Carson + Redução de Kron (ver Objetivo 1)
% =========================================================================
Zabc_mi = [ 0.4576+1j*1.0780,  0.1560+1j*0.5017,  0.1535+1j*0.3849 ;
            0.1560+1j*0.5017,  0.4666+1j*1.0482,  0.1580+1j*0.4236 ;
            0.1535+1j*0.3849,  0.1580+1j*0.4236,  0.4615+1j*1.0651 ]; % Ω/mi

% =========================================================================
% 3. IMPEDÂNCIAS POR TRECHO  [Ω]
%    Comprimentos: trechos 1-2 e 3-4 = 2000 ft; trechos 2-3 e 4-5 = 2500 ft
% =========================================================================
ft2mi = 1 / 5280;   % fator de conversão pés -> milhas

len_ft = [2000, 2500, 2000, 2500];   % comprimentos [ft]: trechos 1-2,2-3,3-4,4-5
len_mi = len_ft * ft2mi;             % comprimentos [mi]

% Célula com as 4 matrizes de impedância de trecho (3x3 cada)
Zseg = cell(4,1);
for k = 1:4
    Zseg{k} = Zabc_mi * len_mi(k);
end

fprintf('=================================================================\n');
fprintf('  PBL1 - BACKWARD-FORWARD SWEEP - IFES Guarapari 2026/1\n');
fprintf('=================================================================\n\n');

fprintf('--- MATRIZES DE IMPEDÂNCIA POR TRECHO [Ω] ---\n');
nomes_seg = {'1-2','2-3','3-4','4-5'};
for k = 1:4
    fprintf('\n  Trecho %s  (%d ft = %.4f mi):\n', nomes_seg{k}, len_ft(k), len_mi(k));
    for i = 1:3
        fprintf('  ');
        for j = 1:3
            fprintf('(%8.5f + j%8.5f)  ', real(Zseg{k}(i,j)), imag(Zseg{k}(i,j)));
        end
        fprintf('\n');
    end
end

% =========================================================================
% 4. POTÊNCIAS DE CARGA POR FASE [VA]
%    Dados: S trifásico / 3 = potência por fase
%    Nó 2: Industrial  - 1000 kVA, fp=0,90 atraso
%    Nó 3: Comercial   -  900 kVA, fp=0,80 atraso
%    Nó 4: Residencial -  750 kVA, fp=0,80 atraso
%    Nó 5: Residencial -  500 kVA, fp=0,80 atraso
% =========================================================================
S3_kVA = [1000, 900, 750, 500];   % potência trifásica [kVA] nos nós 2,3,4,5
fp_carga = [0.90, 0.80, 0.80, 0.80];

% Potência complexa monofásica por fase (S = P + jQ, indutivo => Q > 0)
S1_load = zeros(4,1);   % [VA], complexo
for k = 1:4
    S3  = S3_kVA(k) * 1e3;         % VA
    fp  = fp_carga(k);
    ang = acos(fp);                 % ângulo em rad
    S1_load(k) = (S3/3) * (fp + 1j*sin(ang));  % S = P + jQ por fase
end

% =========================================================================
% 5. FUNÇÃO BFS (definida como função local ao final do script)
%    Chama: bfs(Zseg, V_sub, S1_load, tipo, tol, max_iter)
%    Retorna: V [5x3], I_seg [4x3], I_load [4x3], iters, converged
% =========================================================================

% =========================================================================
% 6. EXECUÇÃO DAS TRÊS SIMULAÇÕES
% =========================================================================
tipos    = {'Z','I','P'};
nomes_nd = {'Industrial','Comercial','Residencial-4','Residencial-5'};

resultados = struct();

for t = 1:3
    tipo = tipos{t};
    [V, I_seg, I_load, iters, conv] = ...
        bfs_sweep(Zseg, V_sub, S1_load, tipo, 1e-6, 100);

    resultados(t).tipo    = tipo;
    resultados(t).V       = V;          % [5 x 3] complexo (nós x fases)
    resultados(t).I_seg   = I_seg;      % [4 x 3] complexo (trechos x fases)
    resultados(t).I_load  = I_load;     % [4 x 3] complexo (nós-carga x fases)
    resultados(t).iters   = iters;
    resultados(t).conv    = conv;
end

% =========================================================================
% 7. IMPRESSÃO DOS RESULTADOS
% =========================================================================
for t = 1:3
    res  = resultados(t);
    tipo = res.tipo;

    fprintf('\n%s\n', repmat('=',1,65));
    fprintf('  CARGA TIPO %s  |  Convergido: %d  |  Iterações: %d\n', ...
            tipo, res.conv, res.iters);
    fprintf('%s\n', repmat('=',1,65));

    % --- Tensões nos nós [p.u.] ---
    fprintf('\n--- TENSÕES NOS NÓS (p.u. em relação à tensão nominal) ---\n');
    for n = 1:5
        V_pu = abs(res.V(n,:)) / Vbase;
        fprintf('  Nó %d: |Va|=%.5f pu  |Vb|=%.5f pu  |Vc|=%.5f pu\n', ...
                n, V_pu(1), V_pu(2), V_pu(3));
    end

    % --- Correntes nos trechos [A] ---
    fprintf('\n--- CORRENTES NOS TRECHOS (magnitudes) ---\n');
    for k = 1:4
        Iabc = abs(res.I_seg(k,:));
        fprintf('  Trecho %s: |Ia|=%.2f A  |Ib|=%.2f A  |Ic|=%.2f A\n', ...
                nomes_seg{k}, Iabc(1), Iabc(2), Iabc(3));
    end

    % --- Correntes de carga [A] ---
    fprintf('\n--- CORRENTES DE CARGA (magnitudes) ---\n');
    for k = 1:4
        Iabc = abs(res.I_load(k,:));
        fprintf('  Nó %d (%s): |Ia|=%.2f A  |Ib|=%.2f A  |Ic|=%.2f A\n', ...
                k+1, nomes_nd{k}, Iabc(1), Iabc(2), Iabc(3));
    end

    % --- Potências nas cargas ---
    fprintf('\n--- POTÊNCIAS CONSUMIDAS POR CARGA ---\n');
    S_total = 0;
    for k = 1:4
        nd = k + 1;   % nó real (2..5)
        % S trifásico = soma das 3 fases: V * conj(I)
        S3 = sum(res.V(nd,:) .* conj(res.I_load(k,:)));
        fprintf('  Nó %d (%s): P=%.2f kW  Q=%.2f kVAr  |S|=%.2f kVA\n', ...
                nd, nomes_nd{k}, real(S3)/1e3, imag(S3)/1e3, abs(S3)/1e3);
        S_total = S_total + S3;
    end

    % --- Balanço de potência ---
    S_fonte  = sum(res.V(1,:) .* conj(res.I_seg(1,:)));
    S_perdas = S_fonte - S_total;
    fprintf('\n--- BALANÇO DE POTÊNCIA ---\n');
    fprintf('  Potência entregue pela fonte:  P=%.2f kW   Q=%.2f kVAr\n', ...
            real(S_fonte)/1e3, imag(S_fonte)/1e3);
    fprintf('  Potência total das cargas:     P=%.2f kW   Q=%.2f kVAr\n', ...
            real(S_total)/1e3, imag(S_total)/1e3);
    fprintf('  Perdas no alimentador:         P=%.2f kW   Q=%.2f kVAr\n', ...
            real(S_perdas)/1e3, imag(S_perdas)/1e3);

    % Armazena para comparação
    resultados(t).S_fonte  = S_fonte;
    resultados(t).S_perdas = S_perdas;
end

% =========================================================================
% 8. TABELA COMPARATIVA FINAL
% =========================================================================
fprintf('\n\n%s\n', repmat('=',1,65));
fprintf('  COMPARAÇÃO DOS TRÊS MODELOS DE CARGA\n');
fprintf('%s\n', repmat('=',1,65));

fprintf('\n%-30s  %10s  %10s  %10s\n', 'Indicador', 'Tipo Z', 'Tipo I', 'Tipo P');
fprintf('%s\n', repmat('-',1,65));

% Tensão mínima no Nó 5, fase A
fprintf('%-30s', 'V mín. Nó 5, Fase A [pu]');
for t = 1:3
    fprintf('  %10.5f', abs(resultados(t).V(5,1)) / Vbase);
end
fprintf('\n');

% Corrente trecho 1-2, fase A
fprintf('%-30s', 'I trecho 1-2, Fase A [A]');
for t = 1:3
    fprintf('  %10.2f', abs(resultados(t).I_seg(1,1)));
end
fprintf('\n');

% Potência entregue pela fonte
fprintf('%-30s', 'P fonte total [kW]');
for t = 1:3
    fprintf('  %10.2f', real(resultados(t).S_fonte)/1e3);
end
fprintf('\n');

% Perdas totais
fprintf('%-30s', 'Perdas totais [kW]');
for t = 1:3
    fprintf('  %10.2f', real(resultados(t).S_perdas)/1e3);
end
fprintf('\n');

% Iterações
fprintf('%-30s', 'Iterações até convergência');
for t = 1:3
    fprintf('  %10d', resultados(t).iters);
end
fprintf('\n\n');

% =========================================================================
% 9. GRÁFICOS COMPARATIVOS
% =========================================================================
nos = 1:5;
cores = [0.13 0.47 0.71 ; 0.17 0.63 0.17 ; 0.84 0.15 0.16];
marcadores = {'o-','s--','d:'};

figure('Name','Perfil de Tensão - Fase A','NumberTitle','off', ...
       'Position',[100 100 820 480]);
hold on; grid on;
for t = 1:3
    V_pu_A = abs(resultados(t).V(:,1)) / Vbase;
    plot(nos, V_pu_A, marcadores{t}, 'Color', cores(t,:), ...
         'LineWidth', 1.8, 'MarkerSize', 8, ...
         'DisplayName', ['Tipo ' resultados(t).tipo]);
end
yline(0.95, 'r--', 'Limite ANEEL (0,95 pu)', 'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
xlabel('Nó', 'FontSize', 11);
ylabel('Tensão [p.u.]', 'FontSize', 11);
title('Perfil de Tensão — Fase A — Modelos Z, I, P', 'FontSize', 12);
legend('Location','southwest','FontSize',10);
xticks(1:5);
xticklabels({'1\newline(Subest.)','2\newline(Industr.)','3\newline(Comerc.)','4\newline(Resid.4)','5\newline(Resid.5)'});
xlim([0.8, 5.2]);

figure('Name','Correntes nos Trechos - Fase A','NumberTitle','off', ...
       'Position',[130 130 820 480]);
hold on; grid on;
trechos_x = [1.5, 2.5, 3.5, 4.5];   % posição no eixo x
for t = 1:3
    I_A = abs(resultados(t).I_seg(:,1));
    plot(trechos_x, I_A, marcadores{t}, 'Color', cores(t,:), ...
         'LineWidth', 1.8, 'MarkerSize', 8, ...
         'DisplayName', ['Tipo ' resultados(t).tipo]);
end
xlabel('Trecho', 'FontSize', 11);
ylabel('Corrente [A]', 'FontSize', 11);
title('Corrente nos Trechos — Fase A — Modelos Z, I, P', 'FontSize', 12);
legend('Location','northeast','FontSize',10);
xticks(trechos_x);
xticklabels({'1-2','2-3','3-4','4-5'});
xlim([1.0, 5.0]);

figure('Name','Potências nas Cargas','NumberTitle','off', ...
       'Position',[160 160 820 500]);
nos_carga = 2:5;
bar_data_P = zeros(4,3);
bar_data_Q = zeros(4,3);
for t = 1:3
    for k = 1:4
        nd = k+1;
        S3 = sum(resultados(t).V(nd,:) .* conj(resultados(t).I_load(k,:)));
        bar_data_P(k,t) = real(S3)/1e3;
        bar_data_Q(k,t) = imag(S3)/1e3;
    end
end

subplot(1,2,1);
b1 = bar(nos_carga, bar_data_P, 0.7);
for t = 1:3, b1(t).FaceColor = cores(t,:); end
xlabel('Nó','FontSize',10); ylabel('P [kW]','FontSize',10);
title('Potência Ativa nas Cargas','FontSize',11);
legend('Tipo Z','Tipo I','Tipo P','Location','northeast','FontSize',9);
grid on; xticks(nos_carga);

subplot(1,2,2);
b2 = bar(nos_carga, bar_data_Q, 0.7);
for t = 1:3, b2(t).FaceColor = cores(t,:); end
xlabel('Nó','FontSize',10); ylabel('Q [kVAr]','FontSize',10);
title('Potência Reativa nas Cargas','FontSize',11);
legend('Tipo Z','Tipo I','Tipo P','Location','northeast','FontSize',9);
grid on; xticks(nos_carga);

fprintf('Simulação concluída. Verifique as figuras geradas.\n\n');

% =========================================================================
% FUNÇÃO LOCAL: bfs_sweep
% =========================================================================
function [V, I_seg, I_load, iters, converged] = ...
         bfs_sweep(Zseg, V_sub, S1_load, tipo, tol, max_iter)
% BFS_SWEEP  Backward-Forward Sweep para alimentador radial trifásico.
%
%  Entradas:
%    Zseg     - cell(4,1) com matrizes 3x3 de impedância por trecho [Ω]
%    V_sub    - vetor [3x1] com fasores de tensão da subestação [V]
%    S1_load  - vetor [4x1] com potência complexa por fase [VA] (nós 2..5)
%    tipo     - 'Z', 'I' ou 'P'
%    tol      - tolerância de convergência [V]
%    max_iter - número máximo de iterações
%
%  Saídas:
%    V        - matriz [5x3] com fasores de tensão (nós x fases) [V]
%    I_seg    - matriz [4x3] com correntes nos trechos (trechos x fases) [A]
%    I_load   - matriz [4x3] com correntes de carga (nós-carga x fases) [A]
%    iters    - número de iterações realizadas
%    converged- flag lógica de convergência

    % Inicialização: perfil plano (todas as tensões = V_sub)
    V       = repmat(V_sub.', 5, 1);   % [5 x 3]
    V_nom   = V_sub.';                 % [1 x 3] tensão nominal de referência
    I_seg   = zeros(4, 3);
    I_load  = zeros(4, 3);
    converged = false;

    for iter = 1:max_iter
        V_old = V;   % guarda tensões da iteração anterior

        % ── BACKWARD SWEEP ────────────────────────────────────────────────
        % Calcula corrente de carga em cada nó (nós 2..5 -> índices k=1..4)
        for k = 1:4
            nd = k + 1;   % nó real
            for ph = 1:3
                V_ph  = V(nd, ph);
                V0_ph = V_nom(ph);
                S1    = S1_load(k);

                switch tipo
                    case 'P'
                        % Potência constante: I = conj(S / V)
                        I_load(k, ph) = conj(S1 / V_ph);

                    case 'I'
                        % Corrente constante: |I| fixo, ângulo rastreia V
                        I_nom = conj(S1 / V0_ph);
                        I_load(k, ph) = abs(I_nom) * exp(1j * angle(V_ph));

                    case 'Z'
                        % Impedância constante: Y = S* / |V_nom|^2
                        Y = conj(S1) / abs(V0_ph)^2;
                        I_load(k, ph) = Y * V_ph;
                end
            end
        end

        % Soma de correntes no sentido da subestação (varredura reversa)
        I_seg(4,:) = I_load(4,:);                        % trecho 4-5
        I_seg(3,:) = I_load(3,:) + I_seg(4,:);          % trecho 3-4
        I_seg(2,:) = I_load(2,:) + I_seg(3,:);          % trecho 2-3
        I_seg(1,:) = I_load(1,:) + I_seg(2,:);          % trecho 1-2

        % ── FORWARD SWEEP ─────────────────────────────────────────────────
        % V(no_k+1) = V(no_k) - Zseg{k} * I_seg(k)
        for k = 1:4
            V(k+1, :) = V(k, :) - (Zseg{k} * I_seg(k,:).').';
        end

        % ── CRITÉRIO DE CONVERGÊNCIA ──────────────────────────────────────
        dV_max = max(max(abs(V(2:5,:) - V_old(2:5,:))));
        if dV_max < tol
            converged = true;
            iters = iter;
            return;
        end
    end

    iters = max_iter;
    warning('BFS nao convergiu em %d iteracoes (tipo %s).', max_iter, tipo);
end
