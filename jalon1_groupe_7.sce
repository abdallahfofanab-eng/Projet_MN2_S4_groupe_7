// =====================================================================
// Projet MN2 - Jalon 1 : Integration numerique
// Sujet Genie civil : energie cumulee a travers une paroi
//   E = integrale de Phi(t) dt   [Wh/m2],  t en heures
//   Phi(t) = U * (T_ext(t) - T_int)   [W/m2]
// Utilisation : exec("jalon1_integration.sce", -1)
// =====================================================================
clear; clc; mode(0);

// ---------------------------------------------------------------------
// 1. Fonctions d'integration (controle des entrees + compteur d'evaluations)
//    f doit etre vectorisee et de la forme y = f(x, p)
// ---------------------------------------------------------------------
function check_inputs(f, a, b, n, nmin, nom)
    if and(type(f) <> [11 13]) then
        error(nom + " : f doit etre une fonction.");
    end
    if ~(isreal(a) & isreal(b)) | size(a, "*") <> 1 | size(b, "*") <> 1 then
        error(nom + " : a et b doivent etre des reels scalaires.");
    end
    if ~(a < b) then
        error(nom + " : il faut a < b.");
    end
    if size(n, "*") <> 1 | ~isreal(n) | n <> round(n) | n < nmin then
        error(nom + " : n doit etre un entier >= " + string(nmin) + ".");
    end
endfunction

function [I, nev] = trapeze(f, a, b, n, p)
    // Trapezes composites sur n sous-intervalles
    check_inputs(f, a, b, n, 1, "trapeze");
    h = (b - a) / n;
    x = a + h * (0:n);
    y = f(x, p);
    nev = size(x, "*");                       // nombre d'evaluations de f
    I = h * (sum(y) - (y(1) + y($)) / 2);
endfunction

function [I, nev] = simpson(f, a, b, n, p)
    // Simpson composite : n pair
    check_inputs(f, a, b, n, 2, "simpson");
    if modulo(n, 2) <> 0 then
        error("simpson : n doit etre pair.");
    end
    h = (b - a) / n;
    x = a + h * (0:n);
    y = f(x, p);
    nev = size(x, "*");
    I = h / 3 * (y(1) + y($) + 4 * sum(y(2:2:$-1)) + 2 * sum(y(3:2:$-2)));
endfunction

// Versions pour donnees echantillonnees (t, y) a pas constant
function I = trapeze_data(t, y)
    h = t(2) - t(1);
    I = h * (sum(y) - (y(1) + y($)) / 2);
endfunction

function I = simpson_data(t, y)
    n = size(t, "*") - 1;
    if modulo(n, 2) <> 0 then
        error("simpson_data : nombre d''intervalles pair requis.");
    end
    h = t(2) - t(1);
    I = h / 3 * (y(1) + y($) + 4 * sum(y(2:2:$-1)) + 2 * sum(y(3:2:$-2)));
endfunction

// ---------------------------------------------------------------------
// 2. Validation sur un polynome de degre <= 3
//    f(x) = x^3 - 2x^2 + x + 1 sur [0,2] ; exact = 8/3
//    Simpson doit etre exact (a la precision machine), pas les trapezes.
// ---------------------------------------------------------------------
function y = poly3(x, p)
    y = x.^3 - 2 * x.^2 + x + 1;
endfunction

Iex = 8 / 3;
[It, nt] = trapeze(poly3, 0, 2, 8, []);
[Is, ns] = simpson(poly3, 0, 2, 8, []);
mprintf("=== Validation (polynome de degre 3, exact = %.12f) ===\n", Iex);
mprintf("Trapezes n=8 : %.12f  erreur = %.2e  (%d evaluations)\n", It, abs(It - Iex), nt);
mprintf("Simpson  n=8 : %.12f  erreur = %.2e  (%d evaluations)\n\n", Is, abs(Is - Iex), ns);

// ---------------------------------------------------------------------
// 3. Modele : paroi, temperatures, flux
// ---------------------------------------------------------------------
Tint = 26;              // [degC] temperature interieure constante
Tm   = 33;              // [degC] temperature exterieure moyenne
Amp  = 7;               // [K]    amplitude journaliere
t0 = 8; t1 = 18;        // [h]    intervalle d'integration
w  = 2 * %pi / 24;      // [rad/h]

// Couches (exterieur -> interieur) : enduit, parpaing creux, enduit
e   = [0.02 0.15 0.015];    // [m]
lam = [1.15 0.75 0.87];     // [W/(m.K)]
Rse = 0.04; Rsi = 0.13;     // [m2.K/W]
RA  = Rse + sum(e ./ lam) + Rsi;     // cas A : sans isolant
RB  = RA + 0.04 / 0.036;             // cas B : + 4 cm polystyrene
UA = 1 / RA;  UB = 1 / RB;           // [W/(m2.K)]
mprintf("=== Parois ===\nU_A = %.4f W/m2.K   U_B = %.4f W/m2.K\n\n", UA, UB);

function T = Text(t)
    T = 33 + 7 * sin(2 * %pi * (t - 9) / 24);     // [degC], t en heures
endfunction

function y = phi(t, p)
    // p = [U, Tint] ; y en W/m2
    y = p(1) * (Text(t) - p(2));
endfunction

function E = E_exacte(t0, t1, p)
    // Primitive analytique [Wh/m2]
    w = 2 * %pi / 24;
    E = p(1) * ((33 - p(2)) * (t1 - t0) + 7 / w * (cos(w * (t0 - 9)) - cos(w * (t1 - 9))));
endfunction

cas_nom = ["A (sans isolant)", "B (avec 4 cm de polystyrene)"];
cas_p   = list([UA Tint], [UB Tint]);
N = [4 8 16 32 64];
tol = 1e-3;                       // tolerance relative visee pour le choix final

// ---------------------------------------------------------------------
// 4. Tableaux de convergence (n = 4, 8, 16, 32, 64)
//    err = |I_n - E_ref| ; ordre observe p = log2(err_{n/2} / err_n)
//    cout = nombre d'evaluations de Phi (n+1)
// ---------------------------------------------------------------------
errT_all = zeros(2, 5); errS_all = zeros(2, 5);
for c = 1:2
    p = cas_p(c);
    Eref = E_exacte(t0, t1, p);
    mprintf("=== Cas %s : E_ref = %.6f Wh/m2 (= %.4f MJ/m2) ===\n", cas_nom(c), Eref, Eref * 3600 / 1e6);
    mprintf("  n | I_trap [Wh/m2] | err_trap [Wh/m2] | p_T | I_simp [Wh/m2] | err_simp [Wh/m2] | p_S | cout\n");
    for k = 1:5
        n = N(k);
        [IT, nT] = trapeze(phi, t0, t1, n, p);
        [IS, nS] = simpson(phi, t0, t1, n, p);
        errT_all(c, k) = abs(IT - Eref);
        errS_all(c, k) = abs(IS - Eref);
        if k == 1 then
            pT = %nan; pS = %nan;
        else
            pT = log2(errT_all(c, k-1) / errT_all(c, k));
            pS = log2(errS_all(c, k-1) / errS_all(c, k));
        end
        mprintf("%3d | %14.6f | %16.3e | %4.2f | %14.6f | %16.3e | %4.2f | %d\n", ..
                n, IT, errT_all(c, k), pT, IS, errS_all(c, k), pS, nT);
    end
    mprintf("\n");
end

// ---------------------------------------------------------------------
// 5. Deux resolutions temporelles des donnees : dt = 1 h et dt = 15 min
//    (en pratique : mesures ou donnees meteo echantillonnees)
// ---------------------------------------------------------------------
mprintf("=== Deux resolutions des donnees (cas A) ===\n");
p = cas_p(1);
Eref = E_exacte(t0, t1, p);
dts = [1, 0.25];
for j = 1:2
    tt = t0:dts(j):t1;
    yy = phi(tt, p);
    ET = trapeze_data(tt, yy);
    ES = simpson_data(tt, yy);
    mprintf("dt = %5.2f h (%2d intervalles) : trap = %.6f  (err rel %.2e) | simp = %.6f  (err rel %.2e)\n", ..
            dts(j), size(tt, "*") - 1, ET, abs(ET - Eref) / Eref, ES, abs(ES - Eref) / Eref);
end
mprintf("\n");

// ---------------------------------------------------------------------
// 6. Decision precision-cout : plus petit n qui atteint l'erreur relative tol
// ---------------------------------------------------------------------
mprintf("=== Choix (tolerance relative %.0e) ===\n", tol);
for c = 1:2
    Eref = E_exacte(t0, t1, cas_p(c));
    kT = find(errT_all(c, :) / Eref <= tol, 1);
    kS = find(errS_all(c, :) / Eref <= tol, 1);
    if isempty(kT) then sT = "non atteint"; else sT = "n = " + string(N(kT)) + " (" + string(N(kT) + 1) + " eval)"; end
    if isempty(kS) then sS = "non atteint"; else sS = "n = " + string(N(kS)) + " (" + string(N(kS) + 1) + " eval)"; end
    mprintf("Cas %s : trapezes -> %s ; Simpson -> %s\n", cas_nom(c), sT, sS);
end

