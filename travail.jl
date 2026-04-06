# ---
# title: Titre du travail
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Auteur
#      prenom: Premier
#      matricule: XXXXXXXX
#      github: premierAuteur
#    - nom: Auteur
#      prenom: Deuxième
#      matricule: XXXXXXXX
#      github: DeuxiAut
# ---

# # Introduction

# # Présentation du modèle
 
# # Implémentation

using CairoMakie
CairoMakie.activate!(px_per_unit=6.0)
using StatsBase, Random
import UUIDs

Random.seed!(2045)

# ----------------------------
# Types
# ----------------------------
Base.@kwdef mutable struct Landscape
    xmin::Int64 = -50
    xmax::Int64 = 50
    ymin::Int64 = -50
    ymax::Int64 = 50
end

Base.@kwdef mutable struct Agent
    x::Int64 = 0
    y::Int64 = 0
    clock::Int64 = 21
    infectious::Bool = false
    vaccinated::Bool = false
    vax_timer::Int64 = 2
    is_detected::Bool = false
    id::UUIDs.UUID = UUIDs.uuid4()
end

Base.@kwdef struct InfectionEvent
    from::UUIDs.UUID
    to::UUIDs.UUID
    time::Int64
    x::Int64
    y::Int64
end

const Population = Vector{Agent}

# ----------------------------
# Random agent generation
# ----------------------------
Random.rand(::Type{Agent}, L::Landscape) = Agent(
    x = rand(L.xmin:L.xmax),
    y = rand(L.ymin:L.ymax)
)

# ----------------------------
# Movement function
# ----------------------------
function move!(A::Agent, L::Landscape; torus=false)
    A.x += rand(Int64(-1):Int64(1))
    A.y += rand(Int64(-1):Int64(1))
    if torus
        A.y = A.y < L.ymin ? L.ymax : A.y
        A.x = A.x < L.xmin ? L.xmax : A.x
        A.y = A.y > L.ymax ? L.ymin : A.y
        A.x = A.x > L.xmax ? L.xmin : A.x
    else
        A.y = clamp(A.y, L.ymin, L.ymax)
        A.x = clamp(A.x, L.xmin, L.xmax)
    end
    return A
end

# ----------------------------
# Population queries
# ----------------------------
isinfectious(agent::Agent) = agent.infectious
ishealthy(agent::Agent) = !agent.infectious

infectious(pop::Population) = filter(isinfectious, pop)
healthy(pop::Population) = filter(ishealthy, pop)
incell(target::Agent, pop::Population) = filter(ag -> (ag.x == target.x && ag.y == target.y), pop)

# ----------------------------
# Simulation
# ----------------------------
function simulation(; maxlength::Int64=1000, population_size::Int64=3750, intervention=true)
    L = Landscape()
    population = [rand(Agent, L) for _ in 1:population_size]
    rand(population).infectious = true

    budget::Float64 = 21000.0
    cout_vax::Float64 = 17.0
    cout_rat::Float64 = 4.0
    morts_totaux::Int64 = 0
    tick::Int64 = 0

    S = Int64[]
    I = Int64[]
    D = Int64[]
    V = Int64[]
    budget_hist = Float64[]
    events = InfectionEvent[]

    while !isempty(infectious(population)) && tick < maxlength
        tick += 1

        # Movement
        for agent in population
            move!(agent, L; torus=false)
        end

        # Infection
        infectieux_du_jour = Random.shuffle(infectious(population))
        for agent in infectieux_du_jour
            neighbors = healthy(incell(agent, population))
            for v in neighbors
                est_protege = v.vaccinated && v.vax_timer <= 0
                if !est_protege && !v.infectious && rand() <= 0.4
                    v.infectious = true
                    v.clock = 21
                    push!(events, InfectionEvent(
                        from = agent.id,
                        to = v.id,
                        time = tick,
                        x = v.x,
                        y = v.y
                    ))
                end
            end
        end

        # Vaccination timer update
        for agent in population
            if agent.vaccinated && agent.vax_timer > 0
                agent.vax_timer -= 1
            end
        end

        # Survival update
        for agent in population
            est_protege = agent.vaccinated && agent.vax_timer <= 0
            if agent.infectious && !est_protege
                agent.clock -= 1
            end
        end

        # Remove dead agents
        morts_step = count(a -> a.clock <= 0, population)
        morts_totaux += morts_step
        population = filter(a -> a.clock > 0, population)

        # Intervention (vaccination + detection)
        if intervention && morts_totaux > 0 && budget > 0 && !isempty(population)
            n = min(Int64(20), Int64(length(population)))
            cibles = StatsBase.sample(population, n; replace=false)

            for a in cibles
                if budget >= cout_rat
                    budget -= cout_rat
                    a.is_detected = a.infectious && rand() <= 0.95
                end
            end

            a_vacciner = filter(a -> !a.infectious && !a.vaccinated, cibles)

            for a in a_vacciner
                if budget >= cout_vax
                    budget -= cout_vax
                    a.vaccinated = true
                    a.vax_timer = 2
                end
            end
        end

        # Store statistics
        push!(S, Int64(length(healthy(population))))
        push!(I, Int64(length(infectious(population))))
        push!(D, morts_totaux)
        push!(V, Int64(count(a -> a.vaccinated, population)))
        push!(budget_hist, budget)
    end

    return (
        tick = tick,
        S = S,
        I = I,
        D = D,
        V = V,
        budget_hist = budget_hist,
        morts_totaux = morts_totaux,
        budget_restant = budget,
        survivants = Int64(length(population)),
        events = events
    )
end

# ----------------------------
# Replicate simulations
# ----------------------------
function replicate_simulations(nrep::Int64; intervention=true)
    results = [simulation(intervention=intervention) for _ in 1:nrep]
    morts = [r.morts_totaux for r in results]
    survivants = [r.survivants for r in results]
    budgets = [21000.0 - r.budget_restant for r in results]
    durations = [r.tick for r in results]

    return (
        results = results,
        morts = morts,
        survivants = survivants,
        budgets = budgets,
        durations = durations,
        mean_morts = mean(morts),
        std_morts = std(morts),
        mean_survivants = mean(survivants),
        std_survivants = std(survivants),
        mean_budget = mean(budgets),
        std_budget = std(budgets),
        mean_duration = mean(durations),
        std_duration = std(durations)
    )
end

# ----------------------------
# Run simulations
# ----------------------------
resultats_sans = simulation(intervention=false)
resultats_avec = simulation(intervention=true)
events = resultats_avec.events

rep_sans = replicate_simulations(Int64(30), intervention=false)
rep_avec = replicate_simulations(Int64(30), intervention=true)

# ----------------------------
# Print summary
# ----------------------------
println("=== COMPARAISON SIMPLE ===")
println("Sans intervention - morts: ", resultats_sans.morts_totaux)
println("Avec intervention - morts: ", resultats_avec.morts_totaux)
println("Réduction de mortalité: ", resultats_sans.morts_totaux - resultats_avec.morts_totaux)
println("Coût de la campagne: ", round(21000.0 - resultats_avec.budget_restant, digits=2))

println()
println("=== RÉPLICATIONS (30) ===")
println("Sans intervention - morts moyens: ", round(rep_sans.mean_morts, digits=2), " ± ", round(rep_sans.std_morts, digits=2))
println("Avec intervention - morts moyens: ", round(rep_avec.mean_morts, digits=2), " ± ", round(rep_avec.std_morts, digits=2))
println("Sans intervention - survivants moyens: ", round(rep_sans.mean_survivants, digits=2), " ± ", round(rep_sans.std_survivants, digits=2))
println("Avec intervention - survivants moyens: ", round(rep_avec.mean_survivants, digits=2), " ± ", round(rep_avec.std_survivants, digits=2))
println("Budget moyen utilisé: ", round(rep_avec.mean_budget, digits=2), " ± ", round(rep_avec.std_budget, digits=2))
println("Durée moyenne avec intervention: ", round(rep_avec.mean_duration, digits=2), " ± ", round(rep_avec.std_duration, digits=2))

# ----------------------------
# Plotting
# ----------------------------
f1 = Figure(size=(900, 600))
ax1 = Axis(f1[1, 1],
    xlabel="Génération",
    ylabel="Population",
    title="Évolution de l'épidémie avec intervention"
)
stairs!(ax1, 1:length(resultats_avec.S), resultats_avec.S, label="Susceptibles", color=:black)
stairs!(ax1, 1:length(resultats_avec.I), resultats_avec.I, label="Infectieux", color=:red)
stairs!(ax1, 1:length(resultats_avec.V), resultats_avec.V, label="Vaccinés", color=:blue)
stairs!(ax1, 1:length(resultats_avec.D), resultats_avec.D, label="Décès cumulés", color=:orange)
axislegend(ax1)
display(f1)

f2 = Figure(size=(900, 600))
ax2 = Axis(f2[1, 1],
    xlabel="Génération",
    ylabel="Budget restant (\$)",
    title="Utilisation du budget pendant l'intervention"
)
lines!(ax2, 1:length(resultats_avec.budget_hist), resultats_avec.budget_hist, color=:green)
display(f2)

f3 = Figure(size=(900, 600))
ax3 = Axis(f3[1, 1],
    xlabel="Répétition",
    ylabel="Morts totaux",
    title="Comparaison des morts: sans vs avec intervention"
)
scatter!(ax3, 1:length(rep_sans.morts), rep_sans.morts, label="Sans intervention")
scatter!(ax3, 1:length(rep_avec.morts), rep_avec.morts, label="Avec intervention")
axislegend(ax3)
display(f3)

if !isempty(events)
    infxn_by_uuid = countmap(getfield.(events, :from))
    nb_inxfn = countmap(collect(values(infxn_by_uuid)))

    f4 = Figure(size=(900, 600))
    ax4 = Axis(f4[1, 1],
        xlabel="Nombre d'infections",
        ylabel="Nombre d'agents",
        title="Distribution du nombre de cas par individu infectieux"
    )
    xvals = sort(collect(keys(nb_inxfn)))
    yvals = [nb_inxfn[x] for x in xvals]
    scatterlines!(ax4, xvals, yvals, color=:black)
    display(f4)

    t = getfield.(events, :time)
    xs = getfield.(events, :x)
    ys = getfield.(events, :y)

    f5 = Figure(size=(900, 600))
    ax5 = Axis(f5[1, 1],
        aspect=1,
        backgroundcolor=:grey97,
        title="Hotspots des infections"
    )
    hm = scatter!(
        ax5,
        xs,
        ys,
        color=t,
        colormap=:viridis,
        strokecolor=:black,
        strokewidth=1,
        colorrange=(0, max(1, resultats_avec.tick)),
        markersize=8
    )
    Colorbar(f5[1, 2], hm, label="Temps d'infection")
    display(f5)

    f6 = Figure(size=(900, 600))
    ax6 = Axis(f6[1, 1],
        xlabel="Temps",
        ylabel="Position x",
        title="Propagation des infections sur l'axe x"
    )
    scatter!(ax6, t, xs, color=:black)
    display(f6)

    f7 = Figure(size=(900, 600))
    ax7 = Axis(f7[1, 1],
        xlabel="Temps",
        ylabel="Position y",
        title="Propagation des infections sur l'axe y"
    )
    scatter!(ax7, t, ys, color=:black)
    display(f7)
end

# Résultat 
# === COMPARAISON SIMPLE ===
# Sans intervention - morts: 2808
# Avec intervention - morts: 443
# Réduction de mortalité: 2365
# Coût de la campagne: 20999.0

# === RÉPLICATIONS (30) ===
# Sans intervention - morts moyens: 2485.57 ± 996.02
# Avec intervention - morts moyens: 627.4 ± 421.22
# Sans intervention - survivants moyens: 1264.43 ± 996.02
# Avec intervention - survivants moyens: 3122.6 ± 421.22
# Budget moyen utilisé: 19347.63 ± 5366.95
# Durée moyenne avec intervention: 748.8 ± 379.38

# # Discussion
# 1. Types & Structs
# Changement : clock augmenté de 20 → 21
# Ajout : vaccinated, vax_timer, is_detected
# Changement : Lattice ±50 au lieu de ±25

# 2. Population Generation
# Changement : génération avec comprehension

# 3. Movement
# Changement : torus=false par défaut, clamp

# 4. Infection
# Ajout : Protection vaccin
# Changement : reset clock à 21

# 5. Survie / Mortalité
# Changement : décrément seulement infectieux non protégés
# Chnagement : suivi morts_totaux

# 6. Intervention
# Nouveau : détection RAT + vaccination + budget

# 7. Statistiques
# Ajout : suivi V_count et budget_count

# 8. Réplications
# Nouveau : réplications + moyennes/écarts-types

# 9. Plotting
# Ajout : plots V/D, budget, hotspots, propagation x/y

"""
# On peut aussi citer des références dans le document `references.bib`, qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# `@` suivi de la clé de citation. Par exemple: @ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.
