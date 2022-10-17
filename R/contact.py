import pandas as pd # Version 0.22.0
import numpy as np # 1.14.0
from scipy.stats import multinomial # 1.0.0
import pymc3 as pm # 3.4.1
import theano # 1.0.1
import theano.tensor as tt
import networkx as nx # 2.1
from collections import OrderedDict
import scipy.stats as stats # 1.0.0
from scipy.optimize import fsolve # 1.0.0


def coerce_index(factor, lookup):
  """
  Convert factor Series into numbers. Similar to the the rethinking package. 
  Useful for fitting random effects.

  Parameters
  ----------
  factor : Series object
    The Series object that will be converted.
  lookup : Series object
    Will use lookup vector to build lookup table. Otherwise uses factor. Just
    pass factor twice if you want to use this as the lookup

  Returns
  -------
  : Series
    factor converted to unique numeric ids.

  """

  lookuptab = pd.Series(lookup.unique()).reset_index()
  lookuptab.columns = ['idnum', factor.name]
  coerced = factor.reset_index().merge(lookuptab, on=factor.name, how="left").set_index(factor.name)

  return((coerced.idnum, lookuptab))


def fixed_effects_model(X, y):
  """
  Build a fixed effects multinomial model in pymc3

  Parameters
  ----------
  X : array
    Design matrix n x p
  y : array or shared theano variable
    Multinomial response variable, n x 3 

  Returns
  -------
  : pymc3 model object containing fitted model
  """

  # Check if y is a shared variable
  if type(y) == tt.sharedvar.TensorSharedVariable:
    yarray = y.eval()
  else:
    yarray = y

  if type(X) == tt.sharedvar.TensorSharedVariable:
    Xmat = X.eval()
  else:
    Xmat = X

  model = pm.Model()

  with model:

    # Transformed coefficients
    betaD = pm.Normal('betad', mu=0, sd=5, shape=(Xmat.shape[1], 1))
    betaI = pm.Normal('betai', mu=0, sd=5, shape=(Xmat.shape[1], 1))

    # Softmax transformation with no contact as baseline
    etad = tt.exp(tt.dot(X, betaD))
    etai = tt.exp(tt.dot(X, betaI))

    p_d = etad / (1 + etad + etai)
    p_i = etai / (1 + etad + etai)
    p_n = 1 - p_d - p_i
    ps = tt.stack([p_n, p_d, p_i], axis=1).reshape((len(yarray), 3)) # Stack column-wise
    #print(ps.tag.test_value.shape)
    #print(betaD.tag.test_value.shape)
    #print(p_d.tag.test_value.shape)
    #print(p_i.tag.test_value.shape)
    #print(p_n.tag.test_value.shape)

    Y_obs = pm.Multinomial("Y_obs", n=1, p=ps, observed=y)
  
  return(model)

def random_effects_model(X, y, a1, a2, spp, spphetero=True, phisd=0.1):
  """
  Multinomial contact model
  
  Random effect of individual species with each species drawing from its own
  distribution.

  Parameters
  ----------
  X : array
    Design matrix n x p
  y : array of theano.shared variable
    Multinomial response variable, n x 3 
  a1 : array
    n x 1 array of unique individual IDs as numbers for first animal in 
    interaction
  a2 : array
    n x 1 array of unique individual IDs as numbers for second animal in the 
    interaction
  spp : array
    For each unique identifier of an individual, the spp identifier.
  phisd : float
    Between species variability in individual heterogeneity.
  spphetero : bool
    If True, fit model with variable individual-level heterogeneity among 
    species. Otherwise, assume similar individual-level heterogeneity for all 
    species. 

  Returns
  -------
  : pymc3 model object
  
  """

  # Check if y or X is a shared variable
  if type(y) == tt.sharedvar.TensorSharedVariable:
    yarray = y.eval()
  else:
    yarray = y

  if type(X) == tt.sharedvar.TensorSharedVariable:
    Xmat = X.eval()
  else:
    Xmat = X

  model = pm.Model()
  num_ids = len(np.unique(np.r_[a1, a2]))
  num_spp = len(np.unique(spp))

  with model:

    # Baseline variance for species-level random effect
    sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))
    sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))

    if spphetero:
      # Proportional perturbance of variance from baseline for other species
      # Direct contact
      phi_vectd = pm.Lognormal("phi_vectd", mu=0, sd=phisd, shape=(num_spp - 1, 1))
      sigma_vectd = phi_vectd*sigma_adirect
      full_sigmad = tt.vertical_stack(sigma_adirect, sigma_vectd)

      # Indirect contact
      phi_vecti = pm.Lognormal("phi_vecti", mu=0, sd=phisd, shape=(num_spp - 1, 1))
      sigma_vecti = phi_vecti*sigma_aindirect
      full_sigmai = tt.vertical_stack(sigma_aindirect, sigma_vecti)

      # Testing the shape 
      # print(sigma_vectd.tag.test_value.shape)
      # print(sigma_direct.tag.test_value.shape)
      # print(full_sigmad.tag.test_value.shape)

      # print(sigma_vecti.tag.test_value.shape)
      # print(sigma_indirect.tag.test_value.shape)
      # print(full_sigmai.tag.test_value.shape)

      # Random effect for each individual, with variance varying by species
      adirect = pm.Normal('adirect', mu=0, sd=full_sigmad[spp], 
                                           shape=(num_ids, 1))
      aindirect = pm.Normal('aindirect', mu=0, sd=full_sigmai[spp], 
                                           shape=(num_ids, 1))
    else:

      # Random effect doesn't vary by species.
      adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                           shape=(num_ids, 1))
      aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                           shape=(num_ids, 1))

    # Fixed effect coefficients
    betaD = pm.Normal('betad', mu=0, sd=2, shape=(Xmat.shape[1], 1))
    betaI = pm.Normal('betai', mu=0, sd=2, shape=(Xmat.shape[1], 1))
    # print(betad.tag.test_value.shape)

    # Softmax transformation with no contact as baseline
    etad = tt.exp(adirect[a1] + adirect[a2] + tt.dot(X, betaD))
    etai = tt.exp(aindirect[a1] + aindirect[a2] + tt.dot(X, betaI))
    p_d = etad / (1 + etad + etai)
    p_i = etai / (1 + etad + etai)
    p_n = 1 - p_d - p_i
    ps = tt.stack([p_n, p_d, p_i], axis=1).reshape((len(yarray), 3)) # Stack column-wise
    # print(ps.tag.test_value.shape)

    # Fit the Multinomial model with random effects of 
    Y_obs = pm.Multinomial("Y_obs", n=1, p=ps, observed=y)
  
  return(model)

def random_effects_model2(X_d, X_id, y, a1, a2, spp, 
                            spphetero=True, phisd=0.1,
                            individual_hetero="both"):
  """
  Multinomial contact model that allows for different design matrices between
  direct and indirect contact and different effect of individual-level
  heterogeneity. 

  Parameters
  ----------
  X_d : array
    Design matrix n x p for direct contacts
  X_id : array
    Design matrix n x p for indirect contacts
  y : array of theano.shared variable
    Multinomial response variable, n x 3 
  a1 : array
    n x 1 array of unique individual IDs as numbers for first animal in the 
    interaction
  a2 : array
    n x 1 array of unique individual IDs as numbers for second animal in the 
    interaction
  spp : array
    An array that is as long as the the number of unique individuals in the 
    analysis, where each item is a number corresponding to which species 
    identity of that individual. Used to generate the species-level effect on
    individual heterogeneity.
  spphetero : bool
    If True, fit model with variable heterogeneity among species. Otherwise,
    assume similar individual-level heterogeneity for all species. 
  individual_hetero : str
    both, neither, direct, or indirect -> Which type of contact should include
    species-level heterogeneity.

  Returns
  -------
  : pymc3 model object
  
  """

  # Check if y or X is a shared variable
  if type(y) == tt.sharedvar.TensorSharedVariable:
    yarray = y.eval()
  else:
    yarray = y

  if type(X_d) == tt.sharedvar.TensorSharedVariable:
    Xmat_d = X_d.eval()
    Xmat_id = X_id.eval()
  else:
    Xmat_d = X_d
    Xmat_id = X_id

  model = pm.Model()
  num_ids = len(np.unique(np.r_[a1, a2]))
  num_spp = len(np.unique(spp))

  with model:

    # Fixed effect coefficients
    betaD = pm.Normal('betad', mu=0, sd=2, shape=(Xmat_d.shape[1], 1))
    betaI = pm.Normal('betai', mu=0, sd=2, shape=(Xmat_id.shape[1], 1))
    # print(betad.tag.test_value.shape)

    # Where is individual-level heterogeneity incorporated?
    if individual_hetero == "both":

      sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))
      sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))

      if spphetero:
        # Proportional perturbance of variance from baseline for other species
        # Direct contact
        phi_vectd = pm.Lognormal("phi_vectd", mu=0, sd=phisd, shape=(num_spp - 1, 1))
        sigma_vectd = phi_vectd*sigma_adirect
        full_sigmad = tt.vertical_stack(sigma_adirect, sigma_vectd)

        # Indirect contact
        phi_vecti = pm.Lognormal("phi_vecti", mu=0, sd=phisd, shape=(num_spp - 1, 1))
        sigma_vecti = phi_vecti*sigma_aindirect
        full_sigmai = tt.vertical_stack(sigma_aindirect, sigma_vecti)

        # Random effect for each individual, with variance varying by species
        adirect = pm.Normal('adirect', mu=0, sd=full_sigmad[spp], 
                                             shape=(num_ids, 1))
        aindirect = pm.Normal('aindirect', mu=0, sd=full_sigmai[spp], 
                                             shape=(num_ids, 1))
      else:

        # Random effect doesn't vary by species.
        adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                             shape=(num_ids, 1))
        aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                             shape=(num_ids, 1))

      etad = tt.exp(adirect[a1] + adirect[a2] + tt.dot(X_d, betaD))
      etai = tt.exp(aindirect[a1] + aindirect[a2] + tt.dot(X_id, betaI))

    elif individual_hetero == "direct":

      sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))

      if spphetero:
        # Proportional perturbance of variance from baseline for other species
        # Direct contact
        phi_vectd = pm.Lognormal("phi_vectd", mu=0, sd=phisd, shape=(num_spp - 1, 1))
        sigma_vectd = phi_vectd*sigma_adirect
        full_sigmad = tt.vertical_stack(sigma_adirect, sigma_vectd)

        # Random effect for each individual, with variance varying by species
        adirect = pm.Normal('adirect', mu=0, sd=full_sigmad[spp], 
                                             shape=(num_ids, 1))
      else:

        # Random effect doesn't vary by species.
        adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                             shape=(num_ids, 1))
      
      etad = tt.exp(adirect[a1] + adirect[a2] + tt.dot(X_d, betaD))
      etai = tt.exp(tt.dot(X_id, betaI))

    elif individual_hetero == "indirect":

      sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))

      if spphetero:

        # Indirect contact
        phi_vecti = pm.Lognormal("phi_vecti", mu=0, sd=phisd, shape=(num_spp - 1, 1))
        sigma_vecti = phi_vecti*sigma_aindirect
        full_sigmai = tt.vertical_stack(sigma_aindirect, sigma_vecti)

        # Random effect for each individual, with variance varying by species
        aindirect = pm.Normal('aindirect', mu=0, sd=full_sigmai[spp], 
                                             shape=(num_ids, 1))
      else:

        aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                             shape=(num_ids, 1))

      etad = tt.exp(tt.dot(X_d, betaD))
      etai = tt.exp(aindirect[a1] + aindirect[a2] + tt.dot(X_id, betaI))

    else: 

      etad = tt.exp(tt.dot(X_d, betaD))
      etai = tt.exp(tt.dot(X_id, betaI))

    # Softmax transformation with no contact as baseline
    p_d = etad / (1 + etad + etai)
    p_i = etai / (1 + etad + etai)
    p_n = 1 - p_d - p_i
    ps = tt.stack([p_n, p_d, p_i], axis=1).reshape((len(yarray), 3)) # Stack column-wise
    # print(ps.tag.test_value.shape)

    # Fit the Multinomial model with random effects of 
    Y_obs = pm.Multinomial("Y_obs", n=1, p=ps, observed=y)
  
  return(model)


def random_effects_model2_stored_feed(X_d, X_id, X_id2, y, a1, a2):
  """
  Multinomial contact model that allows for different design matrices between
  direct and indirect contact and different effects of individual-level
  heterogeneity. Allows for three types of contact.

  Parameters
  ----------
  X_d : array
    Design matrix n x p for direct contacts
  X_id : array
    Design matrix n x p for indirect contacts of type I
  X_id2 : array
    Design matrix n x p for indirect contacts of type II
  y : array of theano.shared variable
    Multinomial response variable, n x 3 
  a1 : array
    n x 1 array of unique individual IDs as numbers for first animal in 
    interaction
  a2 : array
    n x 1 array of unique individual IDs as numbers for second animal in the 
    interaction
  individual_hetero : str
    both, neither, direct, or indirect -> Which type of contact should include
    species-level heterogeneity.

  Returns
  -------
  : pymc3 model object
  
  """

  # Check if y or X is a shared variable
  if type(y) == tt.sharedvar.TensorSharedVariable:
    yarray = y.eval()
  else:
    yarray = y

  if type(X_d) == tt.sharedvar.TensorSharedVariable:
    Xmat_d = X_d.eval()
    Xmat_id = X_id.eval()
    Xmat_id2 = X_id2.eval()
  else:
    Xmat_d = X_d
    Xmat_id = X_id
    Xmat_id2 = X_id2

  model = pm.Model()
  num_ids = len(np.unique(np.r_[a1, a2]))

  with model:

    # Fixed effect coefficients
    betaD = pm.Normal('betad', mu=0, sd=2, shape=(Xmat_d.shape[1], 1))
    betaI = pm.Normal('betai', mu=0, sd=2, shape=(Xmat_id.shape[1], 1))
    betaI2 = pm.Normal('betai2', mu=0, sd=2, shape=(Xmat_id2.shape[1], 1))
    # print(betad.tag.test_value.shape)

    sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))
    sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))
    #sigma_aindirect2 = pm.HalfNormal('sigma_aindirect2', sd=2, shape=(1, 1))

    # Random effect for each individual, with variance varying by species
    adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                         shape=(num_ids, 1))
    aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                         shape=(num_ids, 1))
    # aindirect2 = pm.Normal('aindirect2', mu=0, sd=sigma_aindirect2, 
    #                                      shape=(num_ids, 1))

    etad = tt.exp(adirect[a1] + adirect[a2] + tt.dot(X_d, betaD))
    etai = tt.exp(aindirect[a1] + aindirect[a2] + tt.dot(X_id, betaI))
    etai2 = tt.exp(tt.dot(X_id2, betaI2))

    # Softmax transformation with no contact as baseline
    p_d = etad / (1 + etad + etai + etai2)
    p_i = etai / (1 + etad + etai + etai2)
    p_i2 = etai2 / (1 + etad + etai + etai2)
    p_n = 1 - p_d - p_i - p_i2
    ps = tt.stack([p_n, p_d, p_i, p_i2], axis=1).reshape((len(yarray), 4)) # Stack column-wise
    # print(ps.tag.test_value.shape)

    # Fit the Multinomial model with random effects of 
    Y_obs = pm.Multinomial("Y_obs", n=1, p=ps, observed=y)
  
  return(model)


def poisson_model(X, y, Z1, Z2, offset, num_ids, individual_hetero="both"):
  """
  Build a Poisson representation of the Multinomial model

  Parameters
  ----------
  X : array-like
    The design matrix of the Poisson model
  y : array-like
    The 0/1 response variable for the model
  Z1 : array-like
    Random effect matrix for individual-level random effects for the first 
    species in the interaction.
  Z2 : array-like
    Random effect matrix for individual-level random effects for the second 
    species in the interaction
  offset : array-like
    The log(overlap time) variable to be used as an offset
  num_ids : int
    The number of unique individuals for which to make random effects
  individual_hetero : str
    Either
      "both": individual heterogeneity in both direct and indirect contact is
              included
      "direct": Only individual heterogeneity in direct contact is included
      "indirect": Only individual heterogeneity in indirect contact is included
      "neither": No individual heterogeneity is included in the model

  Returns
  -------
  : pymc3 Model


  """

  pois_model = pm.Model()

  with pois_model:

    beta = pm.Normal("beta", mu=0, sd=2, shape=(X.shape[1], 1))

    if individual_hetero == "both":

      sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))
      sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))
      
      adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                         shape=(num_ids, 1))
      aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                         shape=(num_ids, 1))
      
      full_a = tt.vertical_stack(adirect, aindirect)
      #print(full_a.tag.test_value.shape)
      

    elif individual_hetero == "indirect":

      sigma_aindirect = pm.HalfNormal('sigma_aindirect', sd=2, shape=(1, 1))
      
      aindirect = pm.Normal('aindirect', mu=0, sd=sigma_aindirect, 
                                         shape=(num_ids, 1))
      
      full_a = aindirect 
      #print(full_a.tag.test_value.shape)

    elif individual_hetero == "direct":

      sigma_adirect = pm.HalfNormal('sigma_adirect', sd=2, shape=(1, 1))
      
      adirect = pm.Normal('adirect', mu=0, sd=sigma_adirect, 
                                         shape=(num_ids, 1))
      
      full_a = adirect

    else:
      pass

    if individual_hetero != "neither":
      logmu = offset + tt.dot(X, beta) + tt.dot(Z1, full_a) + tt.dot(Z2, full_a)
    else:
      logmu = offset + tt.dot(X, beta)

    Yobs = pm.Poisson("Yobs", mu=tt.exp(logmu), observed=y)

  return(pois_model)


def obs_pred_from_sim(datdt, model, trace, samples=200, 
                        network='both', fenceval=0):
  """
  Get observed, predicted, and random networks by farm. This function is 
  specific for the bTB farm analysis.

  Parameters
  ----------
  datdt : DataFrame
    Observed contact dataframe
  model : PyMC3 model object
  trace : PyMC3 trace object from model
  samples : int
    Number of random datasets to predict
  network : str
    'both' - Combine both the indirect and direct network
    'direct' - Just output the direct network
    'indirect' - Just output the indirect network
  """

  animals = ['moo', 'deer', 'raccoon', 'opossum']
  Ysim = pm.sample_ppc(trace, model=model, samples=samples)['Y_obs']

  # Specify the type of network to build
  if network == 'both':

    predcscore = [y[:, 0] == 0 for y in Ysim]
    comp = "(moddat.contact_score > 0)"

  elif network == 'direct':

    predcscore = [y[:, 1] == 1 for y in Ysim]
    comp = "(moddat.contact_score == 1)"

  else:

    predcscore = [y[:, 2] == 1 for y in Ysim]
    comp = "(moddat.contact_score == 2)"

  pred_networks = {}
  for p, farm in enumerate(datdt.farm.unique()):
    
    pred_networks[farm] = []
    
    # Clean the observed network
    fenceind = datdt.fence == fenceval
    redind = (datdt.contact_score > -1) & (datdt.spp1.isin(animals) & datdt.spp2.isin(animals))
    tdat = datdt[redind & fenceind]
    farmind = (tdat.farm == farm) #& (tdat.fence == fenceval)
    moddat = tdat[farmind]
    spppair = np.array(["_".join(np.sort([sp1, sp2])) for sp1, sp2 in zip(moddat.spp1, moddat.spp2)])
    moddat = moddat.assign(spppair=spppair)
    
    # Build observed network
    ov = moddat[eval(comp)]
    Gobs = nx.Graph()
    Gobs.add_nodes_from(np.unique(np.concatenate([moddat.id1, moddat.id2])))
    Gobs.add_edges_from(list(zip(ov.id1, ov.id2)))
    
    # Build simulated network
    pnets = []
    for i in range(len(predcscore)):
        
      predc = moddat[predcscore[i][farmind]]
      Gpred = nx.Graph()
      Gpred.add_nodes_from(np.unique(np.concatenate([moddat.id1, moddat.id2])))
      Gpred.add_edges_from(list(zip(predc.id1, predc.id2)))
      
      pnets.append(Gpred)

    # Build random network
    prands = []
    numnodes = len(Gobs.nodes)
    totedge = (numnodes * (numnodes - 1)) / 2
    numedge = len(Gobs.edges)
    pedge = numedge / totedge
    allnodes = np.unique(np.concatenate([moddat.id1, moddat.id2]))
    
    for i in range(len(predcscore)):
        Grand = nx.fast_gnp_random_graph(numnodes, pedge)
        Grand = nx.relabel_nodes(Grand, {n: lab for n, lab in zip(Grand.nodes, allnodes)})
        prands.append(Grand)

    pred_networks[farm] = (Gobs, pnets, prands)

  return(pred_networks)


def random_rewiring(G):
  """
  Given network G, return a network in which the total number of edges are
  randomized amongst available nodes

  Parameters
  ----------
  G : Networkx graph

  Returns
  -------
  : Gnew
    Randomized network
  """

  nodes = list(G.nodes())
  num_edges = len(G.edges())
  new_edges = [tuple(np.random.choice(nodes, size=2, replace=False)) 
                                              for i in range(num_edges)]
  Gnew = nx.Graph()
  Gnew.add_nodes_from(nodes)
  Gnew.add_edges_from(new_edges)
  return(Gnew)


def network_stats(G, net_stats=['mean_degree', 
                                'assortativity', 
                                'mean_clustering', "mean_betweenness",
                                'max_degree', 'max_betweenness',
                                'transitivity']):
  """
  Given a network G, extract relevant network statistics

  Parameters
  ----------
  G : Network X Graph object
    The network on which to calculated statistics

  Returns
  -------
  : dict
    Dictionary of computed statistics
  """


  # Network distributions
  dd_obs = np.sort([d for n, d in G.degree()])
  cd_obs = [d for n, d in nx.clustering(G).items()]
  bd_obs = [d for n, d in nx.betweenness_centrality(G).items()]

  net_fxns = {"mean_degree": lambda G: np.mean(np.sort([d for n, d in G.degree()])), 
              "assortativity": lambda G: nx.degree_pearson_correlation_coefficient(G), 
              "mean_clustering": lambda G: np.mean([d for n, d in nx.clustering(G).items()]), 
              "mean_betweenness": lambda G: np.mean([d for n, d in nx.betweenness_centrality(G).items()]), 
              "max_degree": lambda G: np.max(np.sort([d for n, d in G.degree()])),
              "max_betweenness": lambda G: np.max([d for n, d in nx.betweenness_centrality(G).items()]),
              "transitivity": lambda G: nx.transitivity(G)}

  obsstats = {stat : net_fxns[stat](G) for stat in net_stats}
  return(obsstats)

def betweenness_by_spp(G):
  """
  Get the average betweenness for each species
  """

  an_lets = {'raccoon': 'R', 'deer': 'D', 
             'cow': 'cow', 'opossum': 'O'}
  a = (pd.Series(nx.betweenness_centrality(G)).reset_index()
        .rename(columns={'index': "id", 0: 'value'}))

  for key, an in an_lets.items():
      a.loc[a.id.str.contains(an), 'id'] = key

  return(a.groupby("id")['value'].mean().reset_index())

def degree_by_spp(G):
  """
  Get the average betweenness for each species
  """

  an_lets = {'raccoon': 'R', 'deer': 'D', 
             'cow': 'cow', 'opossum': 'O'}
  a = (pd.Series(dict(G.degree())).reset_index()
        .rename(columns={'index': "id", 0: 'value'}))

  for key, an in an_lets.items():
      a.loc[a.id.str.contains(an), 'id'] = key

  return(a.groupby("id")['value'].mean().reset_index())


def degree_dist(G):
  """
  Extract the degree distribution from G
  """
  return(np.sort([d for n, d in G.degree()])[::-1])

def betweenness_dist(G):
  """
  Extract the degree distribution from G
  """
  return(np.sort([d for n, d in nx.betweenness_centrality(G).items()])[::-1])

def clustering_dist(G):
  """
  Extract the degree distribution from G
  """
  return(np.sort([d for n, d in nx.clustering(G).items()])[::-1])


## Community-level persistence analysis ##

def sample_contact_matrices(model, fit, cols):
  """
  Sample from the posterior distribution of (model, fit) to generate
  a direct and indirect contact probability matrix
  
  Parameters
  ----------
  model : PyMC3 model object
  fit : PyMC3 trace object
  cols : list of column names
    Columns need to be specific to the Michigan study. deer, moo, opossum,
    and raccoon
  
  Returns
  -------
  : tuple
    (Pairwise direct contact matrix as pd.DataFrame, 
     Pairwise indirect contact matrix as pd.DataFrame)
  """
    
  num = fit['betad'].shape[0]
  rand = np.random.randint(0, num, size=1)[0]
  mats = []
  for nm in ['betad', 'betai']:

    bestparams = pd.DataFrame(fit[nm][rand, :, 0], columns=['vals'],
                                  index=cols)

    animals = ['deer', 'moo', 'opossum', "raccoon"]

    deerdeerfallbev_base = bestparams.filter(regex=".*Intercept.*", axis=0)['vals']
    deerdeerfallbev_base.index = ["C(spppair)[T.deer_deer]"]

    vals = []
    for animal in animals:

      anparams = bestparams.filter(regex=".*{0}.*".format(animal), axis=0)
      anest = (deerdeerfallbev_base.values + 
               anparams.filter(regex="^((?!{0}).)*$".format("farm"), axis=0)['vals'])
      if animal == "deer":
          anest = pd.concat([deerdeerfallbev_base, anest])
      #print(np.exp(anest))
      vals.append(np.exp(anest).values)

    contactmat = pd.DataFrame(np.array(vals), index=animals, columns=animals)
    mats.append(contactmat)
  
  # Convert to probabilities
  probd = mats[0] / (1 + mats[0] + mats[1])
  probi = mats[1] / (1 + mats[0] + mats[1])
  contact_probs = (probd, probi)
  
  return(contact_probs)

def sample_contact_matrices_poisson(fit, cols, samps=1):
  """
  Compute the contact rate matrices from the Poisson model

  Parameters
  ----------
  fit : pymc3 Trace object
    Fit from the full Poisson model
  cols : array-like
    The columns from the fixed effect design matrix
  samps : int
    The number of transition matrices to sample from the posterior

  Returns
  -------
  : list of tuples
    Each list item contains a tuple of transition matrices: 
            no contact, direct contact, indirect contact

  """
  
  # Species pairs
  unq_spppair = np.array(['deer_deer', 'deer_moo', 'deer_opossum', 'deer_raccoon', 
                          'moo_moo', 'moo_opossum', 'moo_raccoon', 
                          'opossum_opossum', 'opossum_raccoon','raccoon_raccoon'])

  beta_names = pm.summary(fit).index[0:129].values # 129 is the number of parameters
  name_map = {xname : bnm.split("_")[2] for bnm, xname in zip(beta_names, cols)}

  param_names =  pd.Series(list(name_map.keys()))

  ind_nocontact = ~param_names.str.contains("direct_contact")
  ind_dcontact = param_names.str.contains("direct_contact") & ~(param_names.str.contains("indirect_contact"))
  ind_idcontact = param_names.str.contains("indirect_contact")

  # First calculate no contact since all rates are relative to that
  base_nm = param_names[0] # The intercept, which represents deer_deer no contact in fall in bev

  # Focus on bev, fall
  nc_base = []
  for spp in unq_spppair:
      
    if spp == "deer_deer":
      param_str = base_nm
    else:
      param_str = "+".join([base_nm, "C(spppair)[T.{0}]".format(spp)])
    
    nc_base.append(param_str)

  # Build the appropriate covariate strings for calculating the rates
  dc_base = []
  ic_base = []
  for i, spp in enumerate(unq_spppair):
      
    if spp == "deer_deer":
        param_str_dc = nc_base[i] + "+contact_fact[T.direct_contact]"
        param_str_ic = nc_base[i] + "+contact_fact[T.indirect_contact]"
    
    else:
      param_str_dc = nc_base[i] + "+contact_fact[T.direct_contact]+" +\
                                   "contact_fact[T.direct_contact]:C(spppair)[T.{0}]".format(spp)
      param_str_ic = nc_base[i] + "+contact_fact[T.indirect_contact]+" +\
                                   "contact_fact[T.indirect_contact]:C(spppair)[T.{0}]".format(spp)
    dc_base.append(param_str_dc)
    ic_base.append(param_str_ic)

  # Calculate direct and indirect rates
  indirect_rates = np.array([np.exp(np.array([fit['beta'][:, np.int(name_map[nm]), :].ravel() 
                                      for nm in c_b.split("+")]).sum(axis=0)) 
                              for c_b in ic_base]).T

  direct_rates = np.array([np.exp(np.array([fit['beta'][:, np.int(name_map[nm]), :].ravel() 
                                      for nm in c_b.split("+")]).sum(axis=0)) 
                              for c_b in dc_base]).T

  nocontact_rates = np.array([np.exp(np.array([fit['beta'][:, np.int(name_map[nm]), :].ravel() 
                                      for nm in c_b.split("+")]).sum(axis=0)) 
                              for c_b in nc_base]).T

  # Build candidate transition matrices
  ans = ['deer', 'moo', 'opossum', 'raccoon']
  all_trans_mats = []
  for samp in range(samps):

    index = np.random.randint(0, indirect_rates.shape[0], size=1)
    samp_d = direct_rates[index, :]
    samp_id = indirect_rates[index, :]
    samp_nc = nocontact_rates[index, :]

    trans_mats = [convert_to_trans(x) for x in [samp_nc, samp_d, samp_id]]
    trans_mats = [pd.DataFrame(x, index=ans, columns=ans) for x in trans_mats]
    all_trans_mats.append(tuple(trans_mats))

  return(all_trans_mats)


def convert_to_trans(x):
  """
  Helper function formatting transition matrix
  """
  trans_mat = np.zeros((4, 4))

  # Indices of lower triangle in column order
  tri_ind = (np.array([0,1,2,3,1,2,3,2,3,3]), np.array([0,0,0,0,1,1,1,2,2,3]))

  trans_mat[tri_ind] = x
  upper = np.triu(trans_mat.T, k=1)
  trans_mat = trans_mat + upper

  return(trans_mat)



def community_r0(cmat, abund, prev):
  """
  Species-level and community-level R0

  Parameters
  ----------
  T : 3 x 3 np.array
    Contact matrix 
    diagonal (i.e. intra-specific transmission)
  abund : 1 x 3 np.array
    The relative abundance/density of each species. The values don't matter, only
    the ratios.
  prev : 1 x 3 np.array
    Equilibrium prevalence for each species. 


  """

  # Calculate species-specific R0 values

  C = (cmat.T / np.diag(cmat)).T
  H = (np.tile(abund, 3).reshape(3, 3).T / abund).T # Relative abundance matrix
  P = (np.tile(prev, 3).reshape(3, 3).T / prev).T # Relative prevalence matrix
  W = C * H
  K = W * P

  # Calculate species-level R0 values
  R0s = 1 / (K.sum(axis=1) * (1 - prev))
  pR = np.prod(R0s)

  # Re-calculate W for community matrix
  C = cmat / np.diag(cmat)
  H = np.tile(abund, 3).reshape(3, 3).T / abund
  W = C * H

  # Calculate community-level R0
  p3 = -1
  p2 = np.sum(R0s)
  p1 = (R0s[0]*R0s[1] - R0s[0]*R0s[2] - R0s[1]*R0s[2] + 
        R0s[1]*R0s[2]*W[2, 1]*W[1, 2] + R0s[0]*R0s[1]*W[1, 0]*W[0, 1] + 
        R0s[0]*R0s[2]*W[2, 0]*W[0, 2])
  p0 = (pR*(1  - W[2, 1]*W[1, 2]) -
        pR*(W[0, 1]*W[1, 0] - W[2, 0]*W[1, 2]*W[0, 1]) +
        pR*(W[1, 0]*W[2, 1]*W[0, 2] - W[2, 0]*W[0, 2]))

  # Find community level R0
  R0_comm = np.max(np.roots([p3, p2, p1, p0]))

  return((R0s, np.real(R0_comm)))

def community_r0_direct_indirect(direct, indirect, abund, prev, death, 
                                  shedding, num=20, show_status=True):
  """
  Estimate the direct and indirect contact R0s of 3 species

  Uses back_R0, but tries multiple different starting values if no convergece
  is achieved.

  Parameters
  ----------
  dcmat : 3 x 3 matrix
    The direct contact matrix (can be relative)
  indirect : array-like
    Length 3. The relative species-specific transmission rates from the
    environmental pool. Assuming all species contact the same pool. 
  abund : array-like
    Length 3. The relative abundance of each species
  prev : array-like
    Length 3: Prevalence of each species
  life : array-like
    Length 3. Relative disease-independent death rates of each species.
  shedding : array-like
    Length 3. Relative shedding of each species into the environmental pool.

  Returns
  -------
  : append(R0s_direct, R0s_indirect)
    The direct and indirect R0s for each species given both types of contact.


  """

  R0_both, status = individual_r0_direct_indirect(direct, indirect, abund, prev, 
                                                  death, 
                                                  shedding, num=num, 
                                                  show_status=show_status)

  if status != 1:
    print("Warning: Convergence not obtained for R0 values")

  # Calculate community R0
  R = R0_both[:3] # direct R0s
  I = R0_both[3:] # indirect R0s

  # W matrix. Ratios of abundance and direct contact
  C = direct / np.diag(direct)
  H = np.tile(abund, 3).reshape(3, 3).T / abund
  W = C * H

  # V matrix. Ratio of abundance and indirect contact
  vmat = np.tile(indirect, 3).reshape(3, 3).T / indirect
  V = vmat * H

  p3 = -1
  p2 = np.sum(R0_both)

  p1 = (R[0] * R[1] * (-1 + W[0, 1] * W[1, 0]) + 
        R[0] * R[2] * (-1 + W[0, 2] * W[2, 0]) +
        R[1] * R[2] * (-1 + W[1, 2] * W[2, 1]) + 
        R[0] * I[1] * (-1 + W[1, 0] * V[0, 1]) +
        R[0] * I[2] * (-1 + W[2, 0] * V[0, 2]) + 
        R[1] * I[0] * (-1 + W[0, 1] * V[1, 0]) + 
        R[1] * I[2] * (-1 + W[2, 1] * V[1, 2]) + 
        R[2] * I[0] * (-1 + W[0, 2] * V[2, 0]) + 
        R[2] * I[1] * (-1 + W[1, 2] * V[2, 1]) + 
        I[0] * I[1] * (-1 + V[0, 1] * V[1, 0]) + 
        I[0] * I[2] * (-1 + V[0, 2] * V[2, 0]))

  # Yes, there should be a zero in the last line. Not a typo
  p0 = (R[0] * R[1] * R[2] * (1 - W[1, 2]*W[2, 1] - W[0, 1]*W[1, 0] - W[0, 2]*W[2, 0] + W[0, 1]*W[1, 2]*W[2, 0] + W[0, 2]*W[1, 0]*W[2, 1]) +
        R[0] * R[1] * I[2] * (1 - W[2, 1]*V[1, 2] - W[0, 1]*W[1, 0] - W[0, 2]*W[2, 0] + W[0, 1]*W[2, 0]*V[1, 2] + W[1, 0]*W[2, 1]*V[0, 2]) + 
        R[0] * R[2] * I[1] * (1 - W[1, 2]*V[2, 1] - W[1, 0]*V[0, 1] - W[2, 0]*V[0, 2] + W[1, 2]*W[2, 0]*V[1, 0] + W[0, 2]*W[1, 0]*V[2, 1]) +
        R[1] * R[2] * I[0] * (1 - W[1, 2]*W[2, 1] - W[0, 1]*V[1, 0] - W[0, 2]*V[2, 0] + W[0, 1]*W[1, 2]*V[2, 0] + W[0, 2]*W[2, 1]*V[1, 0]) + 
        R[1] * I[0] * I[2] * (1 - W[2, 1]*V[1, 2] - V[0, 2]*V[2, 0] + W[2, 1]*V[0, 2]*V[1, 0]) + 
        R[2] * I[0] * I[1] * (1 - W[1, 2]*V[2, 1] - V[0, 1]*V[1, 0] + W[1, 2]*V[0, 1]*V[2, 0]) +
        R[0] * I[1] * I[2] * (0 - W[1, 0]*V[0, 1] - W[2, 0]*V[0, 2] + W[2, 0]*V[0, 1]*V[1, 2] + W[1, 0]*V[0, 2]*V[2, 1]))


  R0_comm = np.max(np.roots([p3, p2, p1, p0]))
  return((R0_both, np.real(R0_comm)))

def individual_r0_direct_indirect(direct, indirect, abund, prev, death, 
                                  shedding, num=20, show_status=True):
  """
  Estimate the direct and indirect contact R0s of 3 species

  Uses back_R0, but tries multiple different starting values if no convergece
  is achieved.

  Parameters
  ----------
  dcmat : 3 x 3 matrix
    The direct contact matrix (can be relative)
  indirect : array-like
    Length 3. The relative species-specific transmission rates from the
    environmental pool. Assuming all species contact the same pool. 
  abund : array-like
    Length 3. The relative abundance of each species
  prev : array-like
    Length 3: Prevalence of each species
  death : array-like
    Length 3. Relative disease-independent death rates of each species.
  shedding : array-like
    Length 3. Relative shedding of each species into the environmental pool.

  Returns
  -------
  : append(R0s_direct, R0s_indirect)
    The direct and indirect R0s for each species given both types of contact.


  """

  # Starting guess for R0
  init = R0_single_mode(prev, direct, shedding, abund)

  # Try to find R0s 
  R0_both, disp, status, message = fsolve(back_R0, init / 2, 
                                  (prev, direct, indirect, shedding, abund, death),
                                  factor=np.min(init), full_output=True)
  
  # If unable to find root, try different starting points.
  if status != 1:
      
    seq = np.linspace(0.1, np.max(init), num=num)
    
    count = 0
    while (status != 1) and (count < len(seq)):
        
      R0_both, disp, status, message = fsolve(back_R0, init / 2, 
                  (prev, direct, indirect, shedding, abund, death),
                  factor=seq[count], full_output=True)

      count += 1

    if show_status:
      print("Number of different starting points {0}, Converged: {1}".format(count, status == 1))

  return((R0_both, status))

def community_r0_sei(cmat, abund, prev, life_to_incubation):
  """
  Species-level and community-level R0

  Parameters
  ----------
  cmat : 3 x 3 np.array
    Contact matrix 
    diagonal (i.e. intra-specific transmission)
  abund : 1 x 3 np.array
    The relative abundance/density of each species. The values don't matter, only
    the ratios.
  prev : 1 x 3 np.array
    Equilibrium prevalence for each species. 
  life_to_incubation : 1 x 3 array
    Ratio of  (1 / average life span) / (1 / average incubation) for each species 
  """

  C = (cmat.T / np.diag(cmat)).T
  H = (np.tile(abund, 3).reshape(3, 3).T / abund).T # Relative abundance matrix
  P = (np.tile(prev, 3).reshape(3, 3).T / prev).T # Relative prevalence matrix
  W = C * H
  K = W * P

  # Calculate species-level R0 values for SEI
  R0s = 1 / (K.sum(axis=1) * (1 - prev - life_to_incubation*prev))
  return(R0s)


def prev_eq(x, beta_dc, l, H, R0_D, R0_I):
  """
  Function for estimating equilibrium prevalence given direct and indirect contact
  
  Parameters
  ----------
  x : array-like
      Prevalences. Length n
  beta_dc : n x n array
      Direct transmission matrix. row i, col j is interpreted as transmission from j to i.
      Only need to be relative transmission values
  l : array-like
      length n. Relative shedding rates of species
  H : array-like
      length n. Relative species density
  R0_D : array-like
      length n. Species-specific R0 estimates for direct contact
  R0_I : array-like
      length n. Species-specific R0 estimates for indirect contact.
  
  Notes
  -----
  Use with fsolve
      
  """
    
  W = (beta_dc.T / beta_dc.diagonal()).T
  
  eqs = []
  for i in range(len(x)):
      
    p = x[i]
    nu = x / p
    w = W[i, :]
    wic = l / l[i]
    eta = H / H[i]
    eq = (1 - p)*R0_D[i]*np.sum(w*nu*eta) + (1 - p)*R0_I[i]*np.sum(wic*nu*eta) - 1
    eqs.append(eq)
      
  return(eqs)


def back_R0(R0s, prev, beta_dc, beta_ic, l, H, mus, fixI=False):
  """
  Function to back calculate direct and indirect R0 with ratio data
  
  Parameters
  ----------
  prev : array-like
      Prevalences. Length n
  beta_dc : n x n array
      Direct transmission matrix. row i, col j is interpreted as transmission from j to i.
      Only need to be relative transmission values
  beta_ic : array-like
      length n. Relative contact rates with environmental pool.
  l : array-like
      length n. Relative shedding rates of species
  H : array-like
      length n. Relative species density
  mus : array-like
      Length n. Relative disease-independent death rates of species.
  
  Notes
  -----
  Use with fsolve. Only set up for two or three
  
  """
  
  # Find the maximum possible R0 values given the system
  upper_bound = R0_single_mode(prev, beta_dc, l, H)

  split = np.int((len(R0s) / 2))

  R0s[R0s <= 0] = 1e-10
  frac = 1# - 1e-5
  R0s[:3] = np.where(R0s[:3] > upper_bound[:3], upper_bound[:3]*frac, R0s[:3]) 
  R0s[3:] = np.where(R0s[3:] > upper_bound[3:], upper_bound[3:]*frac, R0s[3:])
  #print(R0s)

  R0s_D = R0s[:split]
  R0s_I = R0s[split:]

  # Compute relative beta
  W = (beta_dc.T / beta_dc.diagonal()).T
  
  eqs = []
  
  # Calculate prevalence equations
  for i in range(len(R0s_D)):
      
    p = prev[i]
    nu = prev / p
    w = W[i, :]
    wic = l / l[i]
    eta = H / H[i]
    eq = (1 - p)*R0s_D[i]*np.sum(w*nu*eta) + (1 - p)*R0s_I[i]*np.sum(wic*nu*eta) - 1
    eqs.append(eq)
  
  # Calculate R0_D ratios
  for i in range(len(R0s_D) - 1):
    
    # Let R0 go to zero if it wants to
    r0ratio = (R0s_D[i] / R0s_D[i + 1])
    param_ratio = (beta_dc[i, i] / beta_dc[i + 1, i + 1]) * (H[i] / H[i + 1]) * (mus[i + 1] / mus[i])
    eq = (param_ratio - r0ratio)
    eqs.append(eq)
  
  # Calculate R0_I ratios
  for i in range(len(R0s_I) - 1): 
    
    r0ratio = (R0s_I[i] / R0s_I[i + 1])
    param_ratio = (beta_ic[i] / beta_ic[i + 1]) * (H[i] / H[i + 1]) * (mus[i + 1] / mus[i]) * (l[i] / l[i + 1])
    eq = (param_ratio  - r0ratio)
    eqs.append(eq)

  eqs = np.array(eqs)

  # Reduce system to 6 equations
  if len(R0s_D) == 3:

    new_eqs = np.empty(6)
    new_eqs[0] = (eqs[0] + 1) - (eqs[1] + 1)
    new_eqs[1:] = eqs[2:]

  else:
    new_eqs = eqs

  return(np.array(new_eqs))


def R0_single_mode(prev, beta_dc, l, H):
  """ 
  Calculate upper bound R0 estimate with only direct or indirect contact
  
  Parameters
  ----------
  prev : array-like
      Prevalences. Length n
  beta_dc : n x n array
      Direct transmission matrix. row i, col j is interpreted as transmission from j to i.
      Only need to be relative transmission values
  l : array-like
      length n. Relative shedding rates of species
  H : array-like
      length n. Relative species density 
  
  Returns
  -------
  : array
      append(R0's direct, R0's indirect) 
  
  """
    
  W = (beta_dc.T / beta_dc.diagonal()).T
  
  R0s = []
  for i in range(len(prev)):
      
    p = prev[i]
    nu = prev / p
    w = W[i, :]
    wic = l / l[i]
    eta = H / H[i]
    
    R0_D = 1 / ((1 - p)*np.sum(w*nu*eta))
    R0s.append(R0_D)
    
    R0_I = 1 / ((1 - p)*np.sum(wic*nu*eta))
    R0s.append(R0_I)
  
  return(np.append(R0s[::2], R0s[1::2])) # Rearrange so that direct contacts are first

def draw_prev():
  """ Draw deer, opossum, and raccoon prevalence sample """

  prev_deer = (0.04 - 0.018) * np.random.random() + 0.018 # From DMU 452 info
  prev_opossum = (0.10798966 - 0.02463366) * np.random.random() + 0.02463366 # From Walter et al. 2013
  prev_raccoon = (0.09087837 - 0.01311344) * np.random.random() + 0.01311344 # From Walter et al. 2013
  return(np.array([prev_deer, prev_opossum, prev_raccoon]))
        

### From Macroeco ###

def sum_of_squares(obs, pred):
    """
    Sum of squares between observed and predicted data
    Parameters
    ----------
    obs : iterable
        Observed data
    pred : iterable
        Predicted data
    Returns
    -------
    float
        Sum of squares
    Notes
    -----
    The length of observed and predicted data must match.
    """

    return np.sum((np.array(obs) - np.array(pred)) ** 2)


def r_squared(obs, pred, one_to_one=False, log_trans=False):
    """
    R^2 value for a regression of observed and predicted data
    Parameters
    ----------
    obs : iterable
        Observed data
    pred : iterable
        Predicted data
    one_to_one : bool
        If True, calculates the R^2 based on the one-to-one line (see [#]_),
        and if False, calculates the standard R^2 based on a linear regression.
        Default False.
    log_trans : bool
        If True, log transforms obs and pred before R^2 calculation.
    Returns
    -------
    float
        R^2 value
    Notes
    -----
    Using the traditional R^2 to compare the fit of observed and predicted
    values may be misleading as the relationship may not be one-to-one but the
    R^2 value may be quite high. The one-to-one option alleviates this problem.
    Note that with the one-to-one option R^2 can be negative.
    Examples
    --------
    >>> import numpy as np
    >>> import macroeco.compare as comp
    >>> # Generate some data
    >>> x_vals = np.linspace(1, 20, num=100)
    >>> y_vals = np.random.normal(4 + x_vals*2, 1)
    >>> # Standard R^2
    >>> comp.r_squared(x_vals, y_vals)
    0.99336568326291697
    >>> # R^2 about the 1:1 line, will be a poor fit (possibly negative)
    >>> comp.r_squared(x_vals, y_vals, one_to_one=True)
    -6.8621799432144988
    >>> # Generate some other data
    >>> y_vals = np.random.normal(x_vals, 1)
    >>> # Normal R^2
    >>> comp.r_squared(x_vals, y_vals)
    0.97651897660174425
    >>> # R^2 on to the one to one line
    >>> comp.r_squared(x_vals, y_vals, one_to_one=True)
    0.97591430200514639
    References
    ----------
    .. [#]
       White, E., Thibault, K., & Xiao, X. (2012). Characterizing the species
       abundance distributions across taxa and ecosystems using a simple
       maximum entropy model. Ecology, 93(8), 1772-8
    """

    if log_trans:
        obs = np.log(obs)
        pred = np.log(pred)

    if one_to_one:
        r_sq = 1 - (sum_of_squares(obs, pred) /
                    sum_of_squares(obs, np.mean(obs)))
    else:
        b0, b1, r, p_value, se = stats.linregress(obs, pred)
        r_sq = r ** 2

    return r_sq

