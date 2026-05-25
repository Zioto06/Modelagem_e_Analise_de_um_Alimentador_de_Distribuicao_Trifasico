"""
=============================================================================
PBL1 - Distribuição de Energia Elétrica - IFES Guarapari - 2026/1
Grupo 1

Modelagem de Alimentador de Distribuição Trifásico e Fluxo de Potência
Método: Backward-Forward Sweep (BFS)

Cargas:
- Tipo Z: impedância constante
- Tipo I: corrente constante
- Tipo P: potência constante

Correção aplicada:
- Para cargas com fator de potência em atraso, usa-se S = P + jQ.
=============================================================================
"""

import json
from pathlib import Path

import numpy as np


# =============================================================================
# 1. PARÂMETROS DO SISTEMA
# =============================================================================

VLL_kV = 11.4                         # Tensão de linha [kV]
Vbase_kV = VLL_kV / np.sqrt(3)        # Tensão de fase base [kV]
Vbase = Vbase_kV * 1e3                # Tensão de fase base [V]
Scc3_MVA = 150                        # Potência de curto-circuito trifásica [MVA]
freq = 60                             # Frequência [Hz]

# Tensão nominal na subestação, fase-neutro, fasores
# Sequência positiva: Va = 0°, Vb = -120°, Vc = +120°
alpha = np.exp(1j * 2 * np.pi / 3)
V_sub = Vbase * np.array([1.0 + 0j, alpha**2, alpha], dtype=complex)


# =============================================================================
# 2. MATRIZ DE IMPEDÂNCIA DE FASE [Zabc] - POR MILHA
# =============================================================================

Zabc_per_mile = np.array([
    [0.4576 + 1j * 1.0780,  0.1560 + 1j * 0.5017,  0.1535 + 1j * 0.3849],
    [0.1560 + 1j * 0.5017,  0.4666 + 1j * 1.0482,  0.1580 + 1j * 0.4236],
    [0.1535 + 1j * 0.3849,  0.1580 + 1j * 0.4236,  0.4615 + 1j * 1.0651]
], dtype=complex)  # Ω/milha


# =============================================================================
# 3. CONVERSÃO DOS COMPRIMENTOS DOS TRECHOS
# =============================================================================

# Comprimentos em pés
segmentos = {
    '1-2': 2000,
    '2-3': 2500,
    '3-4': 2000,
    '4-5': 2500
}

ft_per_mile = 5280.0

Z_seg = {}
for seg, ft in segmentos.items():
    miles = ft / ft_per_mile
    Z_seg[seg] = Zabc_per_mile * miles


# =============================================================================
# 4. IMPEDÂNCIA DA SUBESTAÇÃO - FONTE EQUIVALENTE
# =============================================================================
# Observação:
# A impedância da fonte é calculada e impressa para referência,
# mas não é inserida no fluxo de potência BFS abaixo.
# O Nó 1 é tratado como barramento de tensão fixa.

XR_ratio = 10.0
Sbase_sub = Scc3_MVA * 1e6
Vbase_sub = VLL_kV * 1e3

Z1_mag = (Vbase_sub ** 2) / Sbase_sub
theta = np.arctan(XR_ratio)

R1_sub = Z1_mag * np.cos(theta)
X1_sub = Z1_mag * np.sin(theta)

Z1_sub = R1_sub + 1j * X1_sub
Z0_sub = 3 * Z1_sub

print("=" * 60)
print("IMPEDÂNCIAS DA FONTE EQUIVALENTE (SUBESTAÇÃO)")
print("=" * 60)
print(f"  |Z1| = {Z1_mag * 1000:.4f} mΩ")
print(f"  Z1 = {R1_sub * 1000:.4f} + j{X1_sub * 1000:.4f} mΩ")
print(f"  Z0 = {Z0_sub.real * 1000:.4f} + j{Z0_sub.imag * 1000:.4f} mΩ")
print()


# =============================================================================
# 5. CARGAS NOS NÓS
# =============================================================================
# Fator de potência em atraso -> carga indutiva -> Q positivo
#
# Nó 2: Industrial  - S = 1000 kVA, fp = 0,90 atraso
# Nó 3: Comercial   - S = 900 kVA,  fp = 0,80 atraso
# Nó 4: Residencial - S = 750 kVA,  fp = 0,80 atraso
# Nó 5: Residencial - S = 500 kVA,  fp = 0,80 atraso

loads_3ph = {
    2: {'S_kVA': 1000, 'fp': 0.90, 'nome': 'Industrial'},
    3: {'S_kVA':  900, 'fp': 0.80, 'nome': 'Comercial'},
    4: {'S_kVA':  750, 'fp': 0.80, 'nome': 'Residencial'},
    5: {'S_kVA':  500, 'fp': 0.80, 'nome': 'Residencial'},
}

# Potência complexa monofásica por fase [VA]
S_load_1ph = {}

for nd, d in loads_3ph.items():
    S3 = d['S_kVA'] * 1e3
    fp = d['fp']
    ang = np.arccos(fp)
    S1 = S3 / 3

    # CORREÇÃO PRINCIPAL:
    # Antes estava: fp - j sen(ang), o que representava carga capacitiva.
    # Para fp em atraso / carga indutiva, o correto é fp + j sen(ang).
    S_load_1ph[nd] = S1 * (fp + 1j * np.sin(ang))


# =============================================================================
# 6. FUNÇÃO BACKWARD-FORWARD SWEEP
# =============================================================================

def bfs(load_type='P', tol=1e-6, max_iter=100):
    """
    Backward-Forward Sweep para alimentador radial trifásico.

    Parâmetros:
        load_type: 'Z', 'I' ou 'P'
        tol: tolerância de convergência [V]
        max_iter: número máximo de iterações

    Retorna:
        V: tensões nos nós
        I_seg: correntes nos trechos
        I_load: correntes de carga
        iteration: número de iterações
        converged: True/False
    """

    nodes = [1, 2, 3, 4, 5]
    load_nodes = [2, 3, 4, 5]

    # Inicialização das tensões com perfil plano
    V = {n: V_sub.copy() for n in nodes}

    # Referência nominal para os modelos Z e I
    V_nom = V_sub.copy()

    def load_current(nd, V_nd):
        """
        Calcula a corrente de carga no nó nd para as três fases.
        """
        S1 = S_load_1ph[nd]
        I_load_node = np.zeros(3, dtype=complex)

        for ph in range(3):
            V_ph = V_nd[ph]
            V0_ph = V_nom[ph]

            if load_type == 'P':
                # Potência constante:
                # S = V * conj(I)
                # I = conj(S / V)
                I_load_node[ph] = np.conj(S1 / V_ph)

            elif load_type == 'I':
                # Corrente constante:
                # Mantém a magnitude nominal da corrente.
                #
                # Esta formulação segue o mesmo critério usado no MATLAB enviado:
                # a magnitude é constante e o ângulo acompanha a tensão.
                # Com isso, o modelo I resulta em Q aproximadamente igual a zero.
                I_nom = np.conj(S1 / V0_ph)
                I_load_node[ph] = np.abs(I_nom) * np.exp(1j * np.angle(V_ph))

            elif load_type == 'Z':
                # Impedância constante:
                # Y = conj(S) / |V_nom|²
                # I = Y * V
                Y = np.conj(S1) / (np.abs(V0_ph) ** 2)
                I_load_node[ph] = Y * V_ph

            else:
                raise ValueError("load_type deve ser 'Z', 'I' ou 'P'.")

        return I_load_node

    converged = False

    for iteration in range(1, max_iter + 1):
        V_old = {n: V[n].copy() for n in nodes}

        # ---------------------------------------------------------------------
        # BACKWARD SWEEP
        # ---------------------------------------------------------------------

        I_load = {nd: load_current(nd, V[nd]) for nd in load_nodes}

        I_seg = {}

        I_seg['4-5'] = I_load[5]
        I_seg['3-4'] = I_load[4] + I_seg['4-5']
        I_seg['2-3'] = I_load[3] + I_seg['3-4']
        I_seg['1-2'] = I_load[2] + I_seg['2-3']

        # ---------------------------------------------------------------------
        # FORWARD SWEEP
        # ---------------------------------------------------------------------

        V[2] = V[1] - Z_seg['1-2'] @ I_seg['1-2']
        V[3] = V[2] - Z_seg['2-3'] @ I_seg['2-3']
        V[4] = V[3] - Z_seg['3-4'] @ I_seg['3-4']
        V[5] = V[4] - Z_seg['4-5'] @ I_seg['4-5']

        # ---------------------------------------------------------------------
        # CRITÉRIO DE CONVERGÊNCIA
        # ---------------------------------------------------------------------

        max_dV = max(np.max(np.abs(V[n] - V_old[n])) for n in load_nodes)

        if max_dV < tol:
            converged = True
            break

    return V, I_seg, I_load, iteration, converged


# =============================================================================
# 7. EXECUÇÃO DAS TRÊS SIMULAÇÕES
# =============================================================================

resultados = {}

for tipo in ['Z', 'I', 'P']:
    V, I_seg, I_load, iters, conv = bfs(load_type=tipo)

    resultados[tipo] = {
        'V': V,
        'I_seg': I_seg,
        'I_load': I_load,
        'iters': iters,
        'converged': conv
    }


# =============================================================================
# 8. FUNÇÃO DE RELATÓRIO
# =============================================================================

def relatorio(tipo, res):
    V = res['V']
    I_seg = res['I_seg']
    I_load = res['I_load']
    iters = res['iters']
    conv = res['converged']

    segs = ['1-2', '2-3', '3-4', '4-5']
    load_nodes = [2, 3, 4, 5]

    nomes = {
        2: 'Industrial',
        3: 'Comercial',
        4: 'Residencial-4',
        5: 'Residencial-5'
    }

    print(f"\n{'=' * 65}")
    print(f"  CARGA TIPO {tipo}  |  Convergido: {conv}  |  Iterações: {iters}")
    print(f"{'=' * 65}")

    # -------------------------------------------------------------------------
    # Tensões
    # -------------------------------------------------------------------------

    print("\n--- TENSÕES NOS NÓS (p.u. em relação à tensão nominal) ---")

    for n in range(1, 6):
        V_pu = np.abs(V[n]) / Vbase

        print(
            f"  Nó {n}: "
            f"|Va|={V_pu[0]:.5f} pu  "
            f"|Vb|={V_pu[1]:.5f} pu  "
            f"|Vc|={V_pu[2]:.5f} pu"
        )

    # -------------------------------------------------------------------------
    # Correntes nos trechos
    # -------------------------------------------------------------------------

    print("\n--- CORRENTES NOS TRECHOS (A, magnitudes) ---")

    for seg in segs:
        Ia, Ib, Ic = I_seg[seg]

        print(
            f"  Trecho {seg}: "
            f"|Ia|={np.abs(Ia):.2f} A  "
            f"|Ib|={np.abs(Ib):.2f} A  "
            f"|Ic|={np.abs(Ic):.2f} A"
        )

    # -------------------------------------------------------------------------
    # Correntes de carga
    # -------------------------------------------------------------------------

    print("\n--- CORRENTES DE CARGA (A, magnitudes) ---")

    for nd in load_nodes:
        Ia, Ib, Ic = I_load[nd]

        print(
            f"  Nó {nd} ({nomes[nd]}): "
            f"|Ia|={np.abs(Ia):.2f} A  "
            f"|Ib|={np.abs(Ib):.2f} A  "
            f"|Ic|={np.abs(Ic):.2f} A"
        )

    # -------------------------------------------------------------------------
    # Potências nas cargas
    # -------------------------------------------------------------------------

    print("\n--- POTÊNCIAS CONSUMIDAS POR CARGA (kW + jkVAr) ---")

    S_total = 0 + 0j
    S_cargas = {}

    for nd in load_nodes:
        S_carga = np.sum(V[nd] * np.conj(I_load[nd]))
        S_cargas[nd] = S_carga
        S_total += S_carga

        print(
            f"  Nó {nd} ({nomes[nd]}): "
            f"P={S_carga.real / 1e3:.2f} kW  "
            f"Q={S_carga.imag / 1e3:.2f} kVAr  "
            f"|S|={abs(S_carga) / 1e3:.2f} kVA"
        )

    # -------------------------------------------------------------------------
    # Balanço de potência
    # -------------------------------------------------------------------------

    S_fonte = np.sum(V[1] * np.conj(I_seg['1-2']))
    S_perdas = S_fonte - S_total

    print("\n--- BALANÇO DE POTÊNCIA ---")
    print(
        f"  Potência entregue pela fonte:  "
        f"P={S_fonte.real / 1e3:.2f} kW  "
        f"Q={S_fonte.imag / 1e3:.2f} kVAr"
    )
    print(
        f"  Potência total das cargas:     "
        f"P={S_total.real / 1e3:.2f} kW  "
        f"Q={S_total.imag / 1e3:.2f} kVAr"
    )
    print(
        f"  Perdas no alimentador:         "
        f"P={S_perdas.real / 1e3:.2f} kW  "
        f"Q={S_perdas.imag / 1e3:.2f} kVAr"
    )

    return {
        'V_pu': {n: np.abs(V[n]) / Vbase for n in range(1, 6)},
        'I_seg_mag': {s: np.abs(I_seg[s]) for s in I_seg},
        'I_load_mag': {nd: np.abs(I_load[nd]) for nd in load_nodes},
        'S_carga': S_cargas,
        'S_total': S_total,
        'S_fonte': S_fonte,
        'S_perdas': S_perdas
    }


# =============================================================================
# 9. IMPRESSÃO DOS RESULTADOS
# =============================================================================

print("\n" + "=" * 65)
print("  PBL1 - FLUXO DE POTÊNCIA - BACKWARD-FORWARD SWEEP")
print("  IFES Campus Guarapari | Distribuição de Energia 2026/1")
print("=" * 65)

print("\n--- MATRIZES DE IMPEDÂNCIA POR TRECHO [Ω] ---")

for seg, ft in segmentos.items():
    miles = ft / ft_per_mile
    print(f"\n  Trecho {seg} ({ft} ft = {miles:.4f} mi):")

    Z = Z_seg[seg]

    for i in range(3):
        row = "  " + "  ".join([
            f"({Z[i, j].real:.5f}+j{Z[i, j].imag:.5f})"
            for j in range(3)
        ])
        print(row)


resumo = {}

for tipo in ['Z', 'I', 'P']:
    resumo[tipo] = relatorio(tipo, resultados[tipo])


# =============================================================================
# 10. COMPARAÇÃO FINAL
# =============================================================================

print("\n\n" + "=" * 65)
print("  COMPARAÇÃO DOS TRÊS MODELOS DE CARGA")
print("=" * 65)

print("\nMenor tensão no Nó 5, fase A [pu]:")

for tipo in ['Z', 'I', 'P']:
    vmin = resumo[tipo]['V_pu'][5][0]
    print(f"  Tipo {tipo}: {vmin:.5f} pu")

print("\nCorrente no trecho 1-2, fase A [A]:")

for tipo in ['Z', 'I', 'P']:
    i12a = resumo[tipo]['I_seg_mag']['1-2'][0]
    print(f"  Tipo {tipo}: {i12a:.2f} A")

print("\nPotência total entregue pela fonte [kW]:")

for tipo in ['Z', 'I', 'P']:
    pf = resumo[tipo]['S_fonte'].real / 1e3
    print(f"  Tipo {tipo}: {pf:.2f} kW")

print("\nPotência reativa total entregue pela fonte [kVAr]:")

for tipo in ['Z', 'I', 'P']:
    qf = resumo[tipo]['S_fonte'].imag / 1e3
    print(f"  Tipo {tipo}: {qf:.2f} kVAr")

print("\nPerdas totais no alimentador [kW]:")

for tipo in ['Z', 'I', 'P']:
    pp = resumo[tipo]['S_perdas'].real / 1e3
    print(f"  Tipo {tipo}: {pp:.2f} kW")

print("\nIterações até a convergência:")

for tipo in ['Z', 'I', 'P']:
    iters = resultados[tipo]['iters']
    print(f"  Tipo {tipo}: {iters}")

print("\nSimulação concluída com sucesso.")


# =============================================================================
# 11. EXPORTAÇÃO DOS RESULTADOS PARA JSON NA PASTA DOWNLOADS
# =============================================================================

dados_export = {}

for tipo in ['Z', 'I', 'P']:
    res = resultados[tipo]

    V = res['V']
    I_s = res['I_seg']
    I_l = res['I_load']

    load_nodes = [2, 3, 4, 5]

    v_pu = {
        str(n): (np.abs(V[n]) / Vbase).tolist()
        for n in range(1, 6)
    }

    i_seg = {
        s: np.abs(I_s[s]).tolist()
        for s in I_s
    }

    i_load = {
        str(nd): np.abs(I_l[nd]).tolist()
        for nd in load_nodes
    }

    s_carga = {}

    for nd in load_nodes:
        S_nd = np.sum(V[nd] * np.conj(I_l[nd]))

        s_carga[str(nd)] = {
            'P_kW': float(S_nd.real / 1e3),
            'Q_kVAr': float(S_nd.imag / 1e3),
            'S_kVA': float(abs(S_nd) / 1e3)
        }

    S_fonte = np.sum(V[1] * np.conj(I_s['1-2']))
    S_total = sum(np.sum(V[nd] * np.conj(I_l[nd])) for nd in load_nodes)
    S_perdas = S_fonte - S_total

    dados_export[tipo] = {
        'V_pu': v_pu,
        'I_seg': i_seg,
        'I_load': i_load,
        'S_carga': s_carga,
        'S_fonte_kW': float(S_fonte.real / 1e3),
        'S_fonte_kVAr': float(S_fonte.imag / 1e3),
        'S_total_kW': float(S_total.real / 1e3),
        'S_total_kVAr': float(S_total.imag / 1e3),
        'S_perdas_kW': float(S_perdas.real / 1e3),
        'S_perdas_kVAr': float(S_perdas.imag / 1e3),
        'iters': int(res['iters']),
        'converged': bool(res['converged'])
    }

caminho_saida = Path.home() / "Downloads" / "resultados_bfs.json"

with open(caminho_saida, "w", encoding="utf-8") as f:
    json.dump(dados_export, f, indent=2, ensure_ascii=False)

print(f"\nDados exportados para: {caminho_saida}")
