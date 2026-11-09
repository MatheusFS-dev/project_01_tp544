"""Reference helpers copied from the colleague's HFSS_Base_Dados.py.

The numerical setup and Au_JC are retained. Geometry adaptation is isolated
in the campaign driver. This module has no import-time solver execution.
"""
import os
import math
import cmath
import time

VERSAO_AEDT = "2025.2"
N_PONTOS, MAX_PASSES, DELTA_S = 401, 15, 0.01
MODOS_AR, MODOS_VIDRO = 4, 10
CASO = "A"
P, H_AU, H_SUB, H_AR = 0.320, 0.030, 1.000, 1.000
C0, EPS0 = 299792458.0, 8.8541878128e-12
AQUI = os.path.dirname(os.path.abspath(__file__))

def construir_caso(hfss, caso, angulo):
    """Build the assignment annulus using the reference's material topology.

    Args:
        hfss: Active HFSS design, with reference materials already installed.
        caso: ``A`` for suspended or ``B`` for the SiO2 half space.
        angulo: Incident angle in degrees, in the xz plane.

    Returns:
        Air body, optional substrate body, top z and bottom z in micrometers.

    Raises:
        ValueError: If the case or angle does not belong to the assignment.
        RuntimeError: If either required geometry subtraction fails.
    """
    if caso not in ("A", "B") or angulo not in (0, 15, 30, 45, 60, 75):
        raise ValueError("unsupported assignment case or angle")
    m = hfss.modeler
    ouro = m.create_cylinder("Z", [0, 0, 0], .090, H_AU, 0, "Au_anelc", "Au_JC")
    furo = m.create_cylinder("Z", [0, 0, 0], .050, H_AU, 0, "Au_anelc_furo", "Au_JC")
    if not ouro or not furo or not m.subtract(ouro, furo, keep_originals=False):
        raise RuntimeError("could not hollow the annulus")
    ouro.solve_inside = True
    z_top = H_AU + H_AR
    z_start = -H_AR if caso == "A" else 0
    ar = m.create_box([-P / 2, -P / 2, z_start], [P, P, z_top - z_start], "Air_Domain", "vacuum")
    if not ar or not m.subtract(ar, [ouro], keep_originals=True):
        raise RuntimeError("could not remove Au volume from the air")
    sub = None
    if caso == "B":
        sub = m.create_box([-P / 2, -P / 2, -H_SUB], [P, P, H_SUB], "SiO2_Domain", "SiO2_reference")
        if not sub:
            raise RuntimeError("could not create the substrate")
    hfss["theta_inc"] = f"{angulo}deg"
    return ar, sub, z_top, -H_AR if caso == "A" else -H_SUB

def log(msg):
    print(msg)
    with open(os.path.join(AQUI, "log_base_dados.txt"), "a", encoding="utf-8") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")

def eps_au(lam_nm):
    e = 1.54 - 1.0 / (143.0**2 * (1.0 / lam_nm**2 + 1j / (14500.0 * lam_nm)))
    for A_, l_, g_ in ((1.27, 470.0, 1900.0), (1.10, 325.0, 1060.0)):
        ph = cmath.exp(-1j * math.pi / 4)
        e += A_ / l_ * (ph / (1 / l_ - 1 / lam_nm - 1j / g_) + ph.conjugate() / (1 / l_ + 1 / lam_nm + 1j / g_))
    return e

def face_em(obj, eixo, valor, tol=1e-4):
    k = "xyz".index(eixo)
    for f in obj.faces:
        c = f.center
        if abs(c[k] - valor) < tol and abs(f.normal[k]) > 0.99:
            return f.id
    raise RuntimeError("face %s = %.3f nao encontrada em %s" % (eixo, valor, obj.name))

def fronteiras_e_portas(hfss, ar, sub, z_top, z_bot):
    for obj in [o for o in (ar, sub) if o is not None]:
        for eixo in "xy":
            hfss.assign_lattice_pair([face_em(obj, eixo, -P / 2), face_em(obj, eixo, +P / 2)],
                                     False, "UseScanAngle", "0deg", "theta_inc")
    portas = []
    nmb = MODOS_VIDRO if CASO == "B" else MODOS_AR
    for nome, obj, z, nm in (("FloquetCima", ar, z_top, MODOS_AR),
                             ("FloquetBaixo", sub if CASO == "B" else ar, z_bot, nmb)):
        fz = face_em(obj, "z", z)
        portas.append(hfss.create_floquet_port(fz, [-P / 2, -P / 2, z], [P / 2, -P / 2, z], [-P / 2, P / 2, z], nm, nome))
    return portas

def setup_e_varredura(hfss):
    mod = hfss.odesign.GetModule("AnalysisSetup")
    mod.InsertSetup("HfssDriven", ["NAME:Setup1", "SolveType:=", "Single", "Frequency:=", "666.2THz",
        "MaxDeltaS:=", DELTA_S, "MaximumPasses:=", MAX_PASSES, "MinimumPasses:=", 1,
        "MinimumConvergedPasses:=", 1, "PercentRefinement:=", 30, "IsEnabled:=", True,
        "BasisOrder:=", 1, "DoLambdaRefine:=", True, "DoMaterialLambda:=", True, "SetLambdaTarget:=", False,
        "Target:=", 0.3333, "UseMaxTetIncrease:=", False, "PortAccuracy:=", 2, "UseABCOnPort:=", False,
        "SetPortMinMaxTri:=", False, "UseDomains:=", False, "UseIterativeSolver:=", False,
        "SaveRadFieldsOnly:=", False, "SaveAnyFields:=", True])
    mod.InsertFrequencySweep("Setup1", ["NAME:Varredura", "IsEnabled:=", True, "RangeType:=", "LinearCount",
        "RangeStart:=", "187.4THz", "RangeEnd:=", "666.2THz", "RangeCount:=", N_PONTOS,
        "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
        "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250, "InterpMinSolns:=", 0, "InterpMinSubranges:=", 1,
        "ExtrapToDC:=", False, "InterpUseS:=", True, "InterpUsePortImped:=", True, "InterpUsePropConst:=", True,
        "UseDerivativeConvergence:=", False, "InterpDerivTolerance:=", 0.2, "UseFullBasis:=", True,
        "EnforcePassivity:=", True, "PassivityErrorTolerance:=", 0.0001])
    hfss.save_project()

def para_thz(x):
    x = [float(v) for v in x]
    m = max(x)
    if m > 1e12: return [v * 1e-12 for v in x]
    if m > 1e3:  return [v * 1e-3 for v in x]
    return x

def magnitude(d, ex):
    if hasattr(d, "get_expression_data"):
        r = d.get_expression_data(ex, formula="magnitude")
        if isinstance(r, tuple) and len(r) == 2:
            return list(r[0]), [float(v) for v in r[1]]
        return list(d.primary_sweep_values), [float(v) for v in r]
    return list(d.primary_sweep_values), [float(v) for v in d.data_magnitude(ex)]

def espectro(hfss, pc, pb, modo, nmc, nmb):
    expr_r = ["S(%s:%d,%s:%d)" % (pc, m, pc, modo) for m in range(1, nmc + 1)]
    expr_t = ["S(%s:%d,%s:%d)" % (pb, m, pc, modo) for m in range(1, nmb + 1)]
    d = None
    for tentativa in range(3):
        try:
            _ = hfss.setup_names
            d = hfss.post.get_solution_data(expressions=expr_r + expr_t, setup_sweep_name="Setup1 : Varredura",
                                            report_category="Modal Solution Data")
            if d:
                break
        except Exception as e:
            log("     leitura dos parametros S falhou (tentativa %d): %s" % (tentativa + 1, e))
        time.sleep(3)
        hfss.save_project()
    if not d:
        raise RuntimeError("nao consegui ler os parametros S de Setup1 : Varredura")
    fr, out = None, {}
    for lista, nome in ((expr_r, "R"), (expr_t, "T")):
        soma = None
        for ex in lista:
            x, mg = magnitude(d, ex)
            if fr is None:
                fr = para_thz(x)
            soma = [b * b for b in mg] if soma is None else [a + b * b for a, b in zip(soma, mg)]
        out[nome] = soma
    lam = [C0 / (f * 1e12) * 1e9 for f in fr]
    return lam, out["T"], out["R"]






def materiais(hfss):
    hfss.modeler.model_units = "um"
    f = [150e12 + i * 10e12 for i in range(56)]
    e = [eps_au(C0 / fi * 1e9) for fi in f]
    for nome, y in (("eps_Au", [x.real for x in e]),
                    ("sig_Au", [2 * math.pi * fi * EPS0 * x.imag for fi, x in zip(f, e)])):
        try:
            hfss.create_dataset(nome, f, y, is_project_dataset=True)
        except Exception:
            pass
    if "au_jc" not in [k.lower() for k in hfss.materials.material_keys]:
        au = hfss.materials.add_material("Au_JC")
        au.permittivity = "pwl($eps_Au, Freq)"
        au.conductivity = "pwl($sig_Au, Freq)"
