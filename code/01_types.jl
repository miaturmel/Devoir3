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
