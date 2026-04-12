# ---
# title: Simulation d'une campagne de vaccination
# repository: miaturmel/Devoir3
# auteurs:
#    - nom: Gomez Saucedo
#      prenom: Carla Danahe
#      matricule: 20341379
#      github: CarlaGomez1
#    - nom: Modibo Koné
#      prenom: Maimouna
#      matricule: 20234378
#      github: mkone
#    - nom: Turmel
#      prenom: Mia
#      matricule: 20277557
#      github: miaturmel
# ---

# # Introduction
# Une épidémie est définie par une augmentation anormale du nombre de cas d’une même 
# infection au sein d’une population, sur une période donnée (@ifrcepidemics2022). 
# Elle est le résultat de la propagation de germes comme de virus, de bactéries ou encore de 
# parasites qui peuvent se transmettre d’un individu à un autre de manière directe, par 
# exemple : par le contact physique ou par les fluides de la personne infectée ou indirecte, 
# par exemple : l’air, des objets contaminés, la nourriture, l’eau ou encore des vecteurs, 
# soit des insectes ou animaux (@ifrcepidemics2022). La façon par laquelle une maladie se
# propage dépend de plusieurs facteurs, comme : la gravité et la persistance de la maladie, 
# son apparition pour la première fois au sein d’une population, l’apparition de vecteurs et 
# leur quantité, l’état de santé de la population, le contact entre les personnes ainsi que le 
# nombre de personnes étant vaccinés (@ifrcepidemics2022). Dans certaines situations, il 
# existe des individus pouvant être asymptomatiques, mais tout de même contagieux qui
# peuvent transmettre la maladie sans le savoir. Cela constitue un problème qui rend la
# détection et le contrôle d’une épidémie plus difficile. Pour limiter la propagation, il existe 
# des interventions de santé publique possibles pouvant être mises en place, comme la 
# vaccination ou le dépistage. Cependant, ces interventions ont des contraintes liées aux
# ressources limitées, ce qui impose de faire des choix stratégiques pour leur 
# utilisation.

# Dans cette simulation, la maladie étudiée est transmissible, asymptomatique et toujours
# fatale en l’absence de protection. Cette maladie touche une population de 3740 individus
# n’ayant jamais été confrontée à la maladie étudiée, ils ne sont donc pas immunisés. Donc,
# dans ce contexte biologique, comment l’utilisation de tests de dépistage et de la 
# vaccination, sous contrainte budgétaire, influence-t-elle la propagation de l’épidémie et 
# le nombre total de décès dans la population atteinte? Cette simulation permettra donc de
# modéliser la propagation spatiale d’une infection ainsi que l’utilité des méthodes
# d’intervention, comme le dépistage et les vaccins. La première hypothèse posée sur ce
# modèle est que, sans aucune intervention, l’épidémie pourrait se propager rapidement à
# l’ensemble de la population, entraînant une mortalité très élevée parce que cette maladie
# est asymptomatique et fatale lorsqu’elle est contractée. La deuxième hypothèse est que
# l’introduction d’interventions (dépistages et vaccination) permettra de réduire
# significativement propagation et mortalité de la population. En effet, le dépistage va
# permettre d’identifier les individus infectieux et la vaccination protégera les individus
# sains contre l’infection. La troisième hypothèse est que l’efficacité des interventions va
# fortement dépendre de la gestion du budget de 21000$ disponible pour les vaccins (17$)
# et le dépistage (4$) pour une population de 3750 individus. En effet, si le budget était
# réparti équitablement entre tous les membres de cette population, chaque personne
# n’aurait que 5,60$, juste assez pour payer un seul dépistage. Dans la simulation
# d’épidémies, on s’attend donc à voir une plus grande mortalité dans la population
# lorsqu’il y a peu ou pas d’interventions et, au contraire, lorsque les interventions
# augmentent, le taux de propagation et mortalité diminue dans la population.


# # Présentation du modèle

# ## Suppositions du modèle :
# -La population initiale n’est pas immunisé
# -Le taux d’infection est de 0.4 
# -La durée de la maladie est de 21 jours et est toujours fatale 
# -Les individus infectieux sont asymptomatiques. 
# -La transmission dépend uniquement des contacts entre les individus de la population fermée 
# -Le vaccin ne s’active qu’après deux jours suivant son inoculation
# -Le vaccin protège complètement un individu de la maladie 
# -Les tests de dépistage antigéniques ont une efficacité de 95%, mais ne permettent pas de savoir depuis quand un individu est infectieux
# -Les ressources sont limitées par un budget fixe de 21000$ pour l’ensemble de la population. 
# -Les vaccins coûtent 17$ chacun 
# -Les tests de dépistage coutent 4$ chacun
# -La seule façon de connaître la prévalence de la maladie est par le moyen des tests 

# ## Notre modèle de vaccination :
# Une intervention dans la population n'est possible que s’il y a au moins 1 mort résultant de la maladie dans la population, 
# s’il reste de l’argent dans le budget pour réaliser des tests et administrer des vaccins et si la population n’est pas vide,
# c'est-à-dire si pas tous agents sont morts de la maladie. On sélectionne ensuite un maximum de 20 agents aléatoirement et sans remise dans
# la population à chaque pas de temps. Ainsi, un même individu ne peut pas être sélectionné plus d'une fois pour un dépistage.
# En effet, à chaque jour, la clinique de vaccination a un nombre limité de rendez-vous disponibles en raison du nombre d'employés limités, 
# de la superficie de la clinique et de la quantité de ressources disponibles, par exemple. Si le budget le permet, nous effectuons ensuite
# un RAT sur ces 20 agents. En filtrant ensuite les agents infectieux et ceux ayant déjà reçu un vaccin, nous obtenons un sous-groupe d'individus 
# à partir des 20 sélectionné au départ qui sont sain et non-vacciné et qui sont ceux a qui nous administrerons le vaccins, si le budget le permet.

# # Implémentation

# ## Les packages nécessaires pour simuler le code.

using CairoMakie ## Pour créer des graphiques.
CairoMakie.activate!(px_per_unit=6.0) ## Activation de CairoMakie comme moteur d’affichage graphique.

# px_per_unit=6.0 augmente la résolution des figures.

using StatsBase ## Fournit des fonctions statistiques (tirage aléatoire avec probabilités).
using Random ## Génère des nombres aléatoires.
import UUIDs ## Permet de générer des ID uniques.

# # Code

Random.seed!(2045) ## Garantit des résultats reproductibles.

# Base.@kwdef mutable struct permet de créer une structure mutable qui peut être initialisée avec des valeurs spécifiques.
# Ici, la taille du Landscape est établie.

Base.@kwdef mutable struct Landscape
    xmin::Int64 = -50
    xmax::Int64 = 50
    ymin::Int64 = -50
    ymax::Int64 = 50
end

# Ici, les caractéristiques de départ de l'agent (individus de la population) sont établies.

Base.@kwdef mutable struct Agent
    ## Position
    x::Int64 = 0 
    y::Int64 = 0
    ## Horloge interne de l'agent
    clock::Int64 = 21
    ## Indique si l’agent est infectieux 
    infectious::Bool = false
    ## Indique si l’agent est vacciné
    vaccinated::Bool = false
    ## Indique quand le vaccin prend effet (pour ce code, il prend effet après 2 jours)
    vax_timer::Int64 = 2
    ## Indique si l’agent a été détecté comme infecté
    is_detected::Bool = false
    ## ID unique attribué à chaque agent
    id::UUIDs.UUID = UUIDs.uuid4()
end

# Ici, les caractéristiques de l'événement de transmission de l’infection sont établies.

Base.@kwdef struct InfectionEvent
    ## ID de l’agent ayant transmis l’infection
    from::UUIDs.UUID
    ## ID de l’agent ayant reçu l’infection
    to::UUIDs.UUID
    ## Moment où la transmission a eu lieu dans la simulation
    time::Int64
    ## Position où la transmission a eu lieu
    x::Int64
    y::Int64
end

# Pour ce code, la population est un vecteur contenant plusieurs agents.

const Population = Vector{Agent}

# Ceci permet de générer un agent aléatoire dans un paysage (Landscape) donné
# et de lui assigner une position aléatoire à l’intérieur des limites du Landscape.

Random.rand(::Type{Agent}, L::Landscape) = Agent(
    x = rand(L.xmin:L.xmax),
    y = rand(L.ymin:L.ymax)
)

# Ceci permet de déplacer un agent (A) dans le Landscape (L)

"""
    move!(A::Agent, L::Landscape; torus=false)

Déplace un agent aléatoirement d'un pas de taille 1 dans les deux dimension x et y. 
Le sens du déplacement est aléatoirement défini entre -1 et 1 sur les deux axes.
Si activé, le keyword de la fonction permet que l'agent réapparait de l'autre coté de la lattice
s'il en dépasse les bordures. 

## Arguments et keyword:
A::Agent : identité de l'agent qui subit le déplacement 
L::Landscape : lattice sur laquelle l'agent se déplace
keyword : permet de définir la lattice comme un environnement toroidal. 
    - Si true : l'environnement est toroidal et si un agent dépasse les limite du landscape, 
        il revient de l'autre coté.
    - Si false : l'agent est contraint aux limites du landscape

## Retour:
La fonction retourne la position de l'agent modifiée.
"""
function move!(A::Agent, L::Landscape; torus=false)
    ## Déplacement limité de l'agent
    A.x += rand(Int64(-1):Int64(1))
    A.y += rand(Int64(-1):Int64(1))
    ## Grâce à la fonction torus, si l'agent atteint la bordure du Landscape, il est renvoyé de l’autre côté
    if torus
        A.y = A.y < L.ymin ? L.ymax : A.y
        A.x = A.x < L.xmin ? L.xmax : A.x
        A.y = A.y > L.ymax ? L.ymin : A.y
        A.x = A.x > L.xmax ? L.xmin : A.x
    ## Sinon, l’agent reste à l’intérieur des limites du Landscape
    else
        A.y = clamp(A.y, L.ymin, L.ymax)
        A.x = clamp(A.x, L.xmin, L.xmax)
    end
    return A
end

# Vérifie l'état de l'agent (infectieux ou en santé)

"""
    isinfectious(agent::Agent)

Indique si un agent est infectieux ou non à l'aide de valeurs Booléennes.

## Arguments 
agent::Agent : identité de l'agent dont on vérifie l'état d'infection

## Retour
La fonction retourne une valeur Booléennes : true si l'agent est infectieux, false si non.
"""
isinfectious(agent::Agent) = agent.infectious

"""
    ishealthy(agent::Agent)

Indique si un agent est sain ou non à l'aide de valeurs Booléennes, en comparant l'état de l'agent au contraire de 'agent.infectious'.

## Arguments 
agent::Agent : identité de l'agent dont on vérifie l'état d'infection

## Retour
La fonction retourne une valeur Booléennes : true si agent.infectious == false, false si agent.infectious == true.
"""
ishealthy(agent::Agent) = !agent.infectious

# Filtre la population pour ne garder que les agents infectieux ou sains

"""
    infectious(pop::Population)

Retourne les agents infectieux d'une population en filtrant la population contenant tous les agents pour ne garder que les 
agents infectieux.

## Arguments 
pop::Population = la population contenant tous les agents (infectieux et sains)

## Retour
La fonction retourne une collection contenant uniquement les agents infectieux de la population.
"""
infectious(pop::Population) = filter(isinfectious, pop)

"""
    healthy(pop::Population)

Retourne les agents sains d'une population en filtrant la population contenant tous les agents pour ne garder que les 
agents sains.

## Arguments 
pop::Population = la population contenant tous les agents (infectieux et sains)

## Retour
La fonction retourne une collection contenant uniquement les agents sains de la population.
"""
healthy(pop::Population) = filter(ishealthy, pop)

# incell(target, pop) analyse les agents de la population
# et renvoie seulement ceux qui sont aux mêmes coordonnées x et y

"""
    incell(target::Agent, pop::Population)

Identifie un agent spécifique, puis vérifie les agents qui occupent la même cellule dans le Landscape que cet agent d'intéret.
La fonction retourne ensuite toutes les positions sur les axes x et y de ces agents.

## Arguments 
target::Agent : un agent spécifique à qui la fonction compare la position avec les autres agents de la population
pop::Population : la population contenant tous les agents (infectieux et sains)

## Retour
La fonction retourne une collection contenant uniquement les agents dont les coordonnées sur les axes x et y sont les même
    que l'agent target.
"""
incell(target::Agent, pop::Population) = filter(ag -> (ag.x == target.x && ag.y == target.y), pop)

# Crée le Landscape où vivent les agents

"""
    simulation(; maxlength::Int64=1000, population_size::Int64=3750, intervention=true)

## Description 
1. La fonction établie le Landscape dans lequel la population évolue.
2. Génération de la population contenant tous les agents au départ. 
3. Échantillonnage d'un individu spécifique qui sera le premier individu infectieux.
4. Établissement du budget de départ et des couts relier aux interventions.
5. À chaque tick, ou pas de temps, les agents de la population subissent un déplacement, 
    des agents infectieux sont choisi, la propagation de la maladie est générée et une
    horloge interne de 21 jours est initié qui décompte le temps avant qu'un agent infectieux
    meurt.
6. Les agents qui meurt à la fin du délai de 21 jours sont supprimés de la population.
7. L'intervention est simulée : une groupe de 20 agents sont sélectionnés aléatoirement,
    ils se font tester et certains se font vacciner.

## Arguments ou keywords 
maxlength::Int64=1000 : nombre maximal de ticks ou pas de temps simulés
population_size::Int64=3750 : taille de la population initiale 
intervention=true : active ou non l'intervention. Si true, l'intervention est activée et les RAT et la vaccination sont mis en place.
    Si false, l'intervention n'a pas lieu.

## Retour
La fonction retourne retourne un NamedTuple qui contient:
    - tick = le nombre total de ticks simulés
    - S = historique desvagents sains
    - I = historique des agents infectieux
    - D = historique des morts cumulées
    - V = historique des agents vaccinés
    - budget_hist = l'historique du budget
    - morts_totaux = nombre total de morts
    - budget_restant = montant restant du budget initial
    - survivants = nombre d'agents restant dans la population
    - events = historique des évènements d'infection
"""
function simulation(; maxlength::Int64=1000, population_size::Int64=3750, intervention=true)
    L = Landscape()
    ## Un agent est choisi au hasard pour qu'il soit infectieux au départ
    population = [rand(Agent, L) for _ in 1:population_size]
    rand(population).infectious = true

    ## Budget et coûts de départ
    budget::Float64 = 21000.0
    cout_vax::Float64 = 17.0
    cout_rat::Float64 = 4.0
    morts_totaux::Int64 = 0
    tick::Int64 = 0

    ## Stocke les infos de la simulation à chaque étape (nombres de sains, infectés, morts, etc. à chaque tick)
    S = Int64[]
    I = Int64[]
    D = Int64[]
    V = Int64[]
    budget_hist = Float64[]
    events = InfectionEvent[]

    ## Boucle principale de la simulation (un tick correspond à une étape)
    while !isempty(infectious(population)) && tick < maxlength
        tick += 1 

        ## Déplacement des agents
        for agent in population
            move!(agent, L; torus=false)
        end

        ## Établir un ordre aléatoire des infectieux
        infectieux_du_jour = Random.shuffle(infectious(population))
        for agent in infectieux_du_jour
            ## Déterminer tous les voisins sains sur la même cellule
            neighbors = healthy(incell(agent, population))
            for v in neighbors
                est_protege = v.vaccinated && v.vax_timer <= 0
                if !est_protege && !v.infectious && rand() <= 0.4
                    ## Si non protégé, non infectieux et probabilité < 0.4, l'agent devient infectieux
                    v.infectious = true
                    ## Durée de l'infection et enregistrement de l'événement
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

        ## Délai avant l'effet du vaccin
        for agent in population
            if agent.vaccinated && agent.vax_timer > 0
                agent.vax_timer -= 1
            end
        end

        ## Compte à rebours jusqu'à la mort
        for agent in population
            est_protege = agent.vaccinated && agent.vax_timer <= 0
            if agent.infectious && !est_protege
                agent.clock -= 1
            end
        end

        ## Supprime les agents morts
        morts_step = count(a -> a.clock <= 0, population)
        morts_totaux += morts_step
        population = filter(a -> a.clock > 0, population)

        ## Intervention (vaccination et détection)
        if intervention && morts_totaux > 0 && budget > 0 && !isempty(population)
            n = min(Int64(20), Int64(length(population)))
            cibles = StatsBase.sample(population, n; replace=false)

            ## Détecte les infectieux avec 95 % de chance (test RAT)
            for a in cibles
                if budget >= cout_rat
                    budget -= cout_rat
                    a.is_detected = a.infectious && rand() <= 0.95
                end
            end

            ## Vaccination des agents sains restants (délai de 2 jours avant que le vaccin devienne actif)
            a_vacciner = filter(a -> !a.infectious && !a.vaccinated, cibles)
            for a in a_vacciner
                if budget >= cout_vax
                    budget -= cout_vax
                    a.vaccinated = true
                    a.vax_timer = 2
                end
            end
        end

        ## Enregistre les données à chaque tick
        push!(S, Int64(length(healthy(population))))
        push!(I, Int64(length(infectious(population))))
        push!(D, morts_totaux)
        push!(V, Int64(count(a -> a.vaccinated, population)))
        push!(budget_hist, budget)
    end

    ## Retourne les résultats de la simulation
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

# Fonction pour répéter la simulation plusieurs fois et extraire les données pour comparaison

"""
    replicate_simulations(nrep::Int64; intervention=true)

La fonction effectue plusieurs répétition de la simulation d'infection et collecte des
statistiques pour l'analyse.

## Arguments 
nrep::Int6 : nombre de répétitions de la simulation à effectuer
intervention=true : active ou non l'intervention. Si true, l'intervention est activée et les RAT et la vaccination sont mis en place.
    Si false, l'intervention n'a pas lieu.

## Retour
La fonction retourne retourne un NamedTuple qui contient:
    - results : résultats complets de chaque simulation
    - morts : nombre total de morts par simulation
    - survivants : nombre de survivants par simulation
    - budgets : budgets dépensés par simulation
    - durations : durées de chaque simulation (en nombre de tick)
    - mean_morts : moyenne du nombre de morts
    - std_morts : écart-type du nombre de morts
    - mean_survivants : moyenne du nombre de survivants
    - std_survivants : écart-type des survivants
    - mean_budget : moyenne du budget dépensé
    - std_budget : écart-type du budget dépensé
    - mean_duration : durée moyenne des simulations
    - std_duration : écart-type de la durée des simulations
"""
function replicate_simulations(nrep::Int64; intervention=true)
    results = [simulation(intervention=intervention) for _ in 1:nrep]
    morts = [r.morts_totaux for r in results]
    survivants = [r.survivants for r in results]
    budgets = [21000.0 - r.budget_restant for r in results]
    durations = [r.tick for r in results]
    ## Moyennes et écarts-types pour comparaison
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

# ## Simulations avec et sans intervention et collecte des événements

resultats_sans = simulation(intervention=false);
resultats_avec = simulation(intervention=true);
events = resultats_avec.events;

rep_sans = replicate_simulations(Int64(30), intervention=false);
rep_avec = replicate_simulations(Int64(30), intervention=true);

# ## Interprétation des résultats comparatifs et statistiques

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

# ## Graphiques pour visualiser les résultats

# Évolution de la population

f1 = Figure(size=(900, 600))
ax1 = Axis(f1[1, 1],
    xlabel="Génération",
    ylabel="Population",
    title="Figure 1 : Évolution de l'épidémie avec intervention"
)
stairs!(ax1, 1:length(resultats_avec.S), resultats_avec.S, label="Susceptibles", color=:black)
stairs!(ax1, 1:length(resultats_avec.I), resultats_avec.I, label="Infectieux", color=:red)
stairs!(ax1, 1:length(resultats_avec.V), resultats_avec.V, label="Vaccinés", color=:blue)
stairs!(ax1, 1:length(resultats_avec.D), resultats_avec.D, label="Décès cumulés", color=:orange)
axislegend(ax1)
current_figure()

# Budget restant au fil du temps

f2 = Figure(size=(900, 600))
ax2 = Axis(f2[1, 1],
    xlabel="Génération",
    ylabel="Budget restant (\$)",
    title="Figure 2 : Utilisation du budget pendant l'intervention"
)
lines!(ax2, 1:length(resultats_avec.budget_hist), resultats_avec.budget_hist, color=:green)
current_figure()

# Comparaison des morts sur toutes les réplications

f3 = Figure(size=(900, 600))
ax3 = Axis(f3[1, 1],
    xlabel="Répétition",
    ylabel="Morts totaux",
    title="Figure 3 : Comparaison des morts: sans vs avec intervention"
)
scatter!(ax3, 1:length(rep_sans.morts), rep_sans.morts, label="Sans intervention")
scatter!(ax3, 1:length(rep_avec.morts), rep_avec.morts, label="Avec intervention")
axislegend(ax3)
current_figure()

# Graphiques basés sur les événements d'infection

if !isempty(events)
    infxn_by_uuid = countmap(getfield.(events, :from))
    nb_inxfn = countmap(collect(values(infxn_by_uuid)))

    ## Positions et temps des infections
    t = getfield.(events, :time)
    xs = getfield.(events, :x)
    ys = getfield.(events, :y) 
end;

# Nombre de cas par individu infectieux

f4 = Figure(size=(900, 600))
ax4 = Axis(f4[1, 1],
    xlabel="Nombre d'infections",
    ylabel="Nombre d'agents",
    title="Figure 4 : Distribution du nombre de cas par individu infectieux"
    )
xvals = sort(collect(keys(nb_inxfn)))
yvals = [nb_inxfn[x] for x in xvals]
scatterlines!(ax4, xvals, yvals, color=:black)
current_figure()

# Hotspots des infections

f5 = Figure(size=(900, 600))
ax5 = Axis(f5[1, 1],
     aspect=1,
    backgroundcolor=:grey97,
    title="Figure 5 : Hotspots des infections"
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
current_figure()

# Propagation des infections sur l’axe x

f6 = Figure(size=(900, 600))
ax6 = Axis(f6[1, 1],
    xlabel="Temps",
    ylabel="Position x",
    title="Figure 6 : Propagation des infections sur l'axe x"
    )
scatter!(ax6, t, xs, color=:black)
current_figure()

# Propagation des infections sur l’axe y

f7 = Figure(size=(900, 600))
ax7 = Axis(f7[1, 1],
    xlabel="Temps",
    ylabel="Position y",
    title="Figure 7 : Propagation des infections sur l'axe y"
    )
scatter!(ax7, t, ys, color=:black)
current_figure()

# # Résultats 

# ## Comparaison simple 
# Ce modèle d'épidémie, sans intervention, provoque la mort de 2808 agents, dans une population qui, 
# originellement, en comptais 3750. Si on active l'intervention dans la population, ce nombre de morts 
# diminue à 443. Cela représente une diminution de 2365 agents qui ont succombés à la maladie. En revanche,
# cette campagne représente une dépense de 20 999$. 

# ## Avec réplications
# Si on ajoute des réplications à la simulation, il est possible de voir l'effet de la stochasticité sur 
# la simulation et ainsi de calculer des statistiques tel que des moyennes et des écarts-types nos résultats.
# Par exemple, avec 30 réplications, la simulation sans intervention provoque la mort de (2485.57 ± 996.02) agents
# en moyenne. Les survivants moyens de cette simulation sans intervention sont ainsi de (1264.43 ± 996.02). Pour la
# simulation avec intervention, après 30 réplications, le nombre de morts moyens est de (627.4 ± 421.22), ce qui nous donne 
# (3122.6 ± 421.22) survivants en moyenne. Le budget mmoyen utilisé pour cette intervention est de (19347.63 ± 5366.95)$ et
# la durée moyenne de l'intervention est de (748.8 ± 379.38) ticks ou pas de temps, avec un maximum de 1000. Cette durée
# s'arrête lorsque tous les agents de la population ont succombés à la maladie.

# ## Analyse des graphiques
# ### Figure 1 : Évolution de l'épidémie avec intervention
# Il est possible de voir qu'au début, tous les agents de la population sont susceptibles, aucun ne sont infectieux, 
# morts ou vaccinés. La situation restent ainsi pour les dix premières générations environ. Ensuite, on voit l’effet de 
# l'intervention prendre forme: Les individus sont vaccinés à un rythme très élevé (une centaine de génération environ) 
# jusqu’à l’atteinte d’un plateau, correspondant à l’épuisement du budget. On en comprend aussi que le début de la campagne 
# de vaccination correspond à la première mort dans la population : la ligne orange dans le graphique connait une croissance 
# plus rapide au début (100 premières générations environ). Cette tendance à la croissance positive est suivie par la courbe 
# d’individus infectieux. Par la suite (après les 100 premieres générations), on voit que les décès cumulés continus d’augmenter 
# à rythme certes plus lent qu'au départ, tandis que les individus infectieux dans la population sont quasiment inexistant.

# ### Figure 2 : Utilisation du budget pendant l'intervention
# La figure 2 nous montre justement le moment ou le budget est épuisé. La première utilisation du budget, qui se situe au moins 
# à la 21e génération correspond, comme on le sait, à la première mort d’un agent suite à la maladie, puisque c'est à ce moment 
# que l'intervention est déclenchée. Si on assume que le premier agent infecté est infecté dès la première génération, alors il 
# serait mort à la 21e génération. Le budget s’épuise alors très rapidement et est à sec vers la 100e génération, ce qui confirme 
# les résultats observés dans la figure 1. 

# ### Figure 3 : Comparaison des morts : sans vs avec intervention
# Cette figure compare le nombre total de mort générée dans la simulation avec et sans intervention. Il est possible d'observer
# que généralement, le nombre de mort est beaucoup plus élevé sans intervention qu’avec intervention : l'écart entre les deux 
# totaux de morts se situe généralement autour de 1000 individus. En revanche, le nombre de mort dans la simulatiion sans intervention
# ne dépasse généralement pas le 3000 individus. Ainsi, la population n’est jamais anéantie totalement. De plus, il est possible 
# d’observer qu’à quatre reprises, la simulation sans intervention n’a généré aucune mort.

# ### Figure 4 : Distribution du nombre de cas par individu infectieux
# Cette figure montre la distribution du nombre de cas initié par chaque individus, c'eat-à dire qu'on vérifie combien de fois un 
# certain agent a infecté un autre agent. Il est possible d'observer dans la figure que la majorité des agents (environ 150 individus) 
# ont transmis la maladie environ 2 fois. Environ 70 agents l’ont transmis 3 fois, environ 40 agents l’ont transmis 3 fois et environ 
# 5 agents l’ont transmis 5 fois. Seulement 1 agent semble avoir transmis la maladie environ 35 fois, ce qui représente le maximum de 
# fois que la maladie a été transmise.

# ### Figure 5 : Hotspots des infections
# La figure 5 représente la propagation spatio-temporelle de l’épidémie. Il est alors possible de visualiser la position de l’infection
# à travers le temps. Il est possible d'observer que les premières infections sont initiées dans le coin haut-droit de la lattice (vers
# les positions x = 50 et y = 40) et se déplacent ensuite vers le bas de la lattice et vers le coin bas-gauche (vers les positions x = 10
# et x = 50 et y = -20 et y = -40). Le mouvement que l'infection suit est un mouvement un peu en diagonal.

# ### Figure 6 et 7 : Propagation des infections sur l'axe x et y 
# Ces figures détaillent les résultats observés dans la figure 5. Il est possible de voir que, au fil du temps, l'infection subit une
# translation sur l'axe des x, mais restent, pendant toute la durée de la simulation, entre les bornes de x = 5 à x = 50 environ. La
# lattice est composée d'un axe des x qui se prolonge jusqu'à la borne x = -50, mais l'infection ne semble pas s'y propager (voir figure 6).
# Cependant, la maladie se propage de façon plus directionnelle sur l'axe des y (Figure 7) : les agents infectieux au début de la simulation 
# se trouvent du côté positif de l'axe des y, tandis que, vers la fin de la simulation, l'infection est retrouvée du côté négatif de l'axe des y.

# # Discussion

# ## Modifications du code du modèle
# ### 1. Types & Structs
# Initialement, le clock, donc le délai entre l'infection d'un agent et la mort de celui-ci, était de 20 jours. Pour cette simulation, ce
# délai est augmenté à 21 jours. Ainsi, après l'infection d'un agent, il reste en vie pendant 21 jours puis meurt s'il n'est pas guérit. 
# En revanche, dans la simulation présenté ici, la guérison est impossible. Une fois qu'un agent est infectieux, celui-ci ne peut pas 
# redevenir sain; il meurt obligatoirement après 21 jours. 

# De plus, dans le code modifié, l'identié de l'agent contient maintenant une option qui indique si l'agent est vacciné, le temps écoulé
# depuis sa vaccination et indique si l'infection possible de l'agent sera détectée suite à un RAT. Ces changements permettront d'identifier
# les agents qui pourront recevoir les vaccins dans l'intervention simulée et simuler l'innoculation de deux jours qui active le vaccin. 
# En effet, le vaccin qu'on instaure dans cette version du modèle, comme décrit dans l'introduction, permet de protéger des agents sain
# d'une possible infection. Une fois administré, le vaccin n'est effectif que deux jours plus tard. Ainsi, les agents vaccinés depuis moins
# de deux jours sont encore susceptibles aux infections. L'objet "is_detected" permet de présenter le RAT comme un test plutôt réaliste qui 
# ne dépiste pas efficacement la maladie 100% du temps. 

# En effet, le test ne peut pas faire de faux positif, c'eat-à-dire de détecter la maladie si l’agent est sain. Tous les agents sains qui
# se font tester dans notre clinique en ressortent un test négatif. Cette situation représente une situation idéale, mais manque un peu de 
# réalisme, puisque, en pratique, les tests rapides pour la détection de maladie peuvent faire des faux-positif. Par exemple, il a été 
# démontré que les RAT pour la détection de la COVID-19 avait un faible de taux de faux-positif, possiblement lié au moment auquel le test
# a été effectué (trop tôt ou trop tard dans l'infection) ou la qualité de la réalisation du test (@rhodes1995persistence). Le RAT dans la simulation peut, 
# par contre, faire des faux négatif, c'est-à-dire identifier un agent comme sain alors qu'il est en fait infecté. On admet 5% de probabilité 
# pour cette situation dans la simulation. Ce paramètre ajoute du réalisme à la simulation, mais nous empêche, en revanche, d’être efficace à 
# 100% dans notre intervention et d'identifier tous les agents infectieux avec certitude.

# Finalement, la lattice originale dans le code avait des bornes sur les deux axes de -25 à 25, tandis que dans la simulation présentée ici, ces
# bornes sont plutôt de -50 à 50. ????


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

# On peut aussi citer des références dans le document 'references.bib', qui doit
# être au format BibTeX. Les références peuvent être citées dans le texte avec
# '@' suivi de la clé de citation. Par exemple: ermentrout1993cellular -- la
# bibliographie sera ajoutée automatiquement à la fin du document.

# Le format de la bibliographie est American Physics Society, et les références
# seront correctement présentées dans ce format. Vous ne devez/pouvez pas éditer
# la bibliographie à la main.
