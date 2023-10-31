## Warm up: binary dynamics

We start with a quick overview of binary dynamics. Specifically the orbital mechanics of a gravitationally bound system with similar masses
    - I.e. not a planet-star system

These systems are described by the gravitational two-body problem, which as I'm sure you all remember, has an analytical solution. This means that we can perfectly predict the orbit with simple mathematical formulae.

These systems are also stable, which means that as long as there are no external (or internal) changes to the bodies, the system will remain in its orbit indefinitely.

When we assume no other forces, the orbit is described as a Kepler orbit around the binary barycenter, and for these systems you only need two parameters to describe the system: the semi-major axis, and the eccentricity. 
    - Note that this is to explain the ORBIT specifically. If you want to also describe the stellar system, you obviously also need the masses (and potentially metallicity and age).

With is in mind, let's move on to a system with three bodies. This system can be described by solving the so-called three-body problem. Now, if you google Three-body problem you'll probably see something like  this. GReat book by the way.

## The three-body problem 

Or, you'll see something like this; a classic chaotic three-body system.

## The three-body problem 

The three-body problem, as opposed to the two-body counterpart has no analytical solutions, which means it can only be solved by numerical integration.
    - This does not mean there are NO (stable) solutions, but rather that these have to be painstakingly found numerically. 
        - Example: figure eight
As already visualized in the previous animation, these systems tend to be chaotic and unstable, which means they are extremely sensitive to initial conditions and they are short-lived
    - Short-lived means that the systems often quickly disintegrate.
This thus begs the question: how can there then be triple stars?

## The hierarcichal configuration {.smaller}

Well, the answer to this is that triple stars are always found in a specific three-body configuration known as a hierarchical triple.

This configuration allows for long periods of stabiliy, which is why these exist in nature.

The hierarchical setup means that we can approximate the system as two nested keplerian orbits: a compact, so-called inner binary consisting of the primary and secondary, and a wide outer orbit consisting of the tertiary and the barycenter of the inner binary.

The introduction of a third body means that the number of parameters required to describe the system also increases. 
Remember that for a binary, we only really need two quantities (plus 2-6 to include stellar)
For a triple, this number shoots to nine for just the orbits, and up to 18 if you include stellar. These include the two orbital separations, eccentricities, arguments of pericenter, angular momentum, and mutual inclination. The mutual inclination is a very important parameters in hier. triples, but can be more difficult to picture intuitively. *EXPLAIN INCLINATION*

Here we see a very simple schematic of a hierarchical triple.
## The hierarcichal configuration

But this doesn't really give the full picture, as it might be quite unintuitive for peole how the orbits actually look. So to get a better idea of this, here are some animations of different hierarchical triple configurations. All of these systems are equal mass, but a different parameter is changed along each axis. Can you tell which ones?

These simulations all had zero mutual inclination, meaning the orbit is planar. To get a better visualization of what the mutual inclination actually means, here are some simulations with different mutual inclinations. If you get motion sickness you might want to look away.


## The von Zeipel-Kozai-Lidov mechanism

- When there are more than two bodies present in a gravitational system, we get something called the von Zeipel-Kozai-Lidov, or Kozai-Lidov, Lidov-Kozai, or whatever it's the name of the day is.
- This mechanism introduces perioc changes in the eccentricity and mutual inclination of the triple. It's essentially an exchange between the two, meaning that the eccentricity periodically gets excited as the inclination gets decreased, and vice-versa.
- This means that completely circular systems with an initial inclination may become highly eccentric. This is likely the most important effect in a triple system, as most triples have a non-zero inclination, meaning that eccentricity excitations are common.
- The KL mechanism is a secular effect, which means the timescales are substantially higher than the dynamic timescales. The effect is stronger the tighter the triple is, or the closer to 90 the inclination is. The excitation is mostly of the inner eccentricity.


- Because of the secular nature of the KL cycle, it is harder to visualize in a video, because you need to simulate longer, which means you need more frames to get a smooth animation, which takes a long time to render. Nevertheless, here is a video showing the inclination between the two orbits changing over larger timescales.

- As I already mentioned, the eccentricity excitation happens mainly in the inner binary, which is hard to see in this video. So in the next one, I plot just the orbit of the two inner bodies with respect to each other.

- We can take an even closer look at the orbital elements themselves as a function of time. Here we see the two eccentricites along with the mutual inclination as a function of time (in number of outer orbits.) As expected we see the periodic changes in the inner eccentricity and mutual inclination, and their inverse relation.

## Stability
Now let's move on to the topic of stability. I have already mentioned that hierarchical triples are usually stable inside a specific set of orbital parameters and mass ratios. However, as you've already seen, triples can and will often destabilize and break the hierarchy. 
- The fact is that any (physical) three-body system is inherently unstable, and will destabilize in a given time, just for some systems this time can be longer than the estimated lifetime of the universe. 
- So how do we know if a triple is stable or not? One of the more popular ways of doing this is using the so-called Mardling-Aarseth stability criterion, as seen here. 
- This is an analytical approximation for the *critical* separation ratio, i.e. the ratio between the outer and the inner semi-major axes where the system will destabilize.
    - If the real sma ratio is equal to or smaller than the critical sma ratio, the system is labeled as unstable.
- This critical ratio can change during the lifetime of the triple due to stellar evolution.
- This means that a triple can go from stable to unstable.
    - For example, if the primary has strong winds, the inner orbit will widen and the semi-major axis ratio will gradually move towards the critical value.
    - Supernova kicks

## What happens to destabilized triples?
So what actually happens when a triple becomes unstable? 
- Well, some of them actually remain hierarchical, and this is simply because the formula presented is, as mentioned, just an approximation, or a fit, to simulated data. It makes straight lines in the parameter space when the real cutoffs are not straight lines.
    - It should be more likely probability
- The ones that do not remain stable might enter what is sometimes called a *democratic* phase, which is essentially what we saw earlier in the animation.
    - This is a phase characterized by chaotic dynamics, and is often short-lived, meaning that it will often quickly disintegrate.
- You might also get extremely close passages due to strong KL effects.
- These scenarios have several potential outcomes
    - One of the stars might get ejected or become a drifter
        - Explain drifter
    - You might get a close passage that results in a collision/merger
    - You can get ionization, in which all three stars disintegrate
    - You can get more exotic phenomena such as a (partial) tidal disruption, if you have a BH and a star in a close passage
- What if left after one of these outcomes?
    - Well you can get a binary in which one of the components is a merger product,
    - you can get one binary and an ejected/drifted single star
    - in the case of ionization you can get three single stars
    - or other exotica

The theme of stability takes us nicely to the topic of my current project, which is about simulating massive triple stars on the edge of stability

## Quick sidestep: secular vs dynamical evolution
- But first, a very quick sidestep to secular vs dynamical evolution
- Does anyone know?
- Essentially these are two ways of simulating triples, where one of the main differences is the timescales.
- For both cases you are solving sets of ODEs, but they are different for the two.
    - For secular evolution, you are solving ODEs for the orbital elements themselves. I.e. a, e, i, ..
        - As you might imagine, these elements change on long timescales, which means you can take large timesteps, which again means you can do very long simulations with extremely cheap computational costs.
- In secular evolution you have no information about the exact positions of the bodies in their orbits, so we say the system is orbit-averaged. 
    - Note that this only works when the triple is in a hierarchical configuration.

- On the other side you have dynamical evolution, in which the positions and velocities of each component is evolved over time. These of course change extremely fast, so you have to take smaller timesteps, and thus the computational cost increases dramatically.


So this leads to my project, where we combine these two techniques to get an overview of the outcomes of massive triples that reach the edge of stability.
## Massive, unstable triples
- For this project, we start by performing secular evolution on a population of massive triples.
    - This was done by Floris using TRES
    - Check out his paper
- These systems are evolved until they reach some stoping condition, one of which is stability.
- I then take these systems and convert them from secular to dynamical systems
    - Remember that secular evolution does not have information about the positions and velocities of the bodies.
    - We only the orbital and stellar parameters at the end of the simulations.
    - So what to do? 
    - Well, when converting to the dynamical code, we randomly place each system in ten uniformly distributed positions in the orbit.
- We can then simulate these systems dynamically until some stopping condition

The summary of outcomes is shown in this figure
## Massive, unstable triples

I also want to quickly talk a bit about the code I have written for doing these dynamical simulations, which is called Syzygy.jl
## Syzygy.jl
  - Julia code for dynamic simulations of hierarchical multiples
  - Supports arbitrary multiplicity and configurations
    - Binary, Triples, 2+2 quad, 3+1 quad, etc... 
  - Very fast: 2x - 60x faster than other nbody codes
    - Also orders of magnitude more stable
  - Uses formalism from Hamers and Portegies Zwart, 2016 as visualized in this figure.

## Syzygy.jl example: triple setup 