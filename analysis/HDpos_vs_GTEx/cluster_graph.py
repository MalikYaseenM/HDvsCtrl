
from itertools import chain, combinations, islice, product
from matplotlib.patches import PathPatch, Circle
from matplotlib.path import Path
from matplotlib.pyplot import *
import networkx as nx
import numpy as np
import pandas

relabels = [
    ('REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS_', 'REACTOME_RESPITORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS...'),
    ('REACTOME_THROMBIN_SIGNALLING_THROUGH_PROTEINASE_ACTIVATED_RECEPTORS_PARS', 'REACTOME_THROMBIN_SIGNALLING...'),
    ('REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING_AND_DOWNSTREAM_TRANSMISSION_IN_THE_POSTSYNAPTIC_CELL', 'REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING...'),
    ('REACTOME_NEF_MEDIATES_DOWN_MODULATION_OF_CELL_SURFACE_RECEPTORS_BY_RECRUITING_THEM_TO_CLATHRIN_ADAPTERS','REACTOME_NEF_MEDIATES_DOWN_MODULATION_OF_CELL_SURFACE_RECEPTORS...'),
    ('REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL','REACTOME_IMMUNOREGULATORY_INTERACTIONS...'),
    ('REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC','REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY...'),
    ('REACTOME_NFKB_AND_MAP_KINASES_ACTIVATION_MEDIATED_BY_TLR4_SIGNALING_REPERTOIRE','REACTOME_NFKB_AND_MAP_KINASES_ACTIVATION...'),
    ('REACTOME_TRAF6_MEDIATED_INDUCTION_OF_NFKB_AND_MAP_KINASES_UPON_TLR7_8_OR_9_ACTIVATION','REACTOME_TRAF6_MEDIATED_INDUCTION_OF_NFKB_AND_MAP_KINASES...'),
    ('REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1','REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20...'),
    ('REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS_THAT_BIND_AU_RICH_ELEMENTS','REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS...'),
    ('REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S','REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX...'),
    ('REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES','REACTOME_METABOLISM_OF_AMINO_ACIDS...'),
    ('REACTOME_PROSTACYCLIN_SIGNALLING_THROUGH_PROSTACYCLIN_RECEPTOR','REACTOME_PROSTACYCLIN_SIGNALING...'),
    ('REACTOME_REGULATION_OF_INSULIN_SECRETION_BY_GLUCAGON_LIKE_PEPTIDE1','REACTOME_REGULATION_OF_INSULIN_SECRETION'),
    ('REACTOME_SEMA4D_INDUCED_CELL_MIGRATION_AND_GROWTH_CONE_COLLAPSE','REACTOME_SEMA4D_INDUCED_CELL_MIGRATION'),
    ('REACTOME_RAS_ACTIVATION_UOPN_CA2_INFUX_THROUGH_NMDA_RECEPTOR','REACTOME_RASACTIVATION_UPON_CA2_INFLUX')
]
def shortened_geneset_names(df) :
    gs_label_map = {_:_ for _ in df.index}
    for k,v in relabels :
        gs_label_map[k] = v
    return [gs_label_map[_] for _ in df.index]

def compute_jaccard(sig_gsea) :
    jacc = pandas.DataFrame(0,columns=sig_gsea.index,index=sig_gsea.index)
    for (k1,_1), (k2,_2) in combinations(sig_gsea.leadingEdge.items(),2) :
        s1, s2 = set(_1.split(' ')), set(_2.split(' '))
        j = 0
        inter = len(s1.intersection(s2))
        if len(s1.union(s2)) != 0 :
            j = inter/len(s1.union(s2))
        jacc.loc[k1,k2] = jacc.loc[k2,k1] = j

    return jacc

def compute_gene_mat(gsea,padj=0.05) :
    print('Significantly enriched: {}'.format(sig_gsea.shape[0]))

    all_leading_edge_genes = set(chain(*[_.split(' ') for _ in sig_gsea.leadingEdge]))

    gene_mat = pandas.DataFrame(0,columns=sig_gsea.index,index=all_leading_edge_genes)
    for k, _ in sig_gsea.leadingEdge.items() :
        genes = _.split(' ')
        gene_mat[k].loc[genes] = 1

    return gene_mat

def get_pos_bounds(pos) :
    min_x, min_y, max_x, max_y = sys.maxsize, sys.maxsize, -sys.maxsize, -sys.maxsize
    for k,p in pos.items() :
        min_x = min(min_x, p[0])
        min_y = min(min_y, p[1])
        max_x = max(max_x, p[0])
        max_y = max(max_y, p[1])
    return min_x, max_x, min_y, max_y

def cluster_graph(jacc,figsize=(30,30),colormap={},ax=None) :

    jacc_g = nx.from_numpy_matrix(jacc.values)
    jacc_g = nx.relabel_nodes(jacc_g,dict(enumerate(jacc.index)))
    
    # 30 nodes per column
    cursor = [0,1e6]

    comps = sorted(nx.connected_component_subgraphs(jacc_g),key=lambda x: len(x.edges))
    comp_pos = []
    comp_bounds = [sys.maxsize, -sys.maxsize, sys.maxsize, -sys.maxsize]
    for subg in comps :
        # linear graphs have only one path
        linear = False
        try :
            next(nx.chain_decomposition(subg))
        except StopIteration :
            linear = True

        if linear : # linear graph, draw vertically
            nodes,_ = zip(*list(subg.nodes.items()))
            pos = {}
            for node in nodes :
                pos[node] = [cursor[0], cursor[1]]
                cursor[1] -= 1
            pos_df = pandas.DataFrame(pos)
        else :
            #pos = nx.spring_layout(subg, weight='weight', seed=1337)
            pos = nx.nx_pydot.graphviz_layout(subg,prog='neato',seed=1337)
            
            pos_df = pandas.DataFrame(pos)

            # get the bounds of the positions
            min_x, min_y = pos_df.min(axis=1)
            max_x, max_y = pos_df.max(axis=1)
            
            # scale the positions so that they are between 0 and factor
            factor = 3
            norm_pos_df = pos_df.sub(pos_df.min(axis=1),axis=0)
            assert np.allclose(norm_pos_df.min(axis=1),[0,0])
            norm_pos_df = factor*norm_pos_df.div(norm_pos_df.max(axis=1),axis=0)
            assert np.allclose(norm_pos_df.max(axis=1),[factor,factor])
            
            # large subgraphs then to be too compact, expand out the size of
            # the graph based on the number of nodes
            if len(subg.nodes) > 9 :
                norm_pos_df *= np.sqrt(len(subg.nodes)-9)
            # 3-mer graphs are a little too big, make them smaller
            elif len(subg.nodes) == 3 :
                norm_pos_df *= 0.5

            # get the bounds of the positions after scaling
            min_x, min_y = norm_pos_df.min(axis=1)
            max_x, max_y = norm_pos_df.max(axis=1)
            
            norm_pos_df = norm_pos_df.add([cursor[0],cursor[1]-max_y],axis=0)
            
            # get the bounds of the positions after scaling
            min_x, min_y = norm_pos_df.min(axis=1)
            max_x, max_y = norm_pos_df.max(axis=1)
            
            assert min_x == cursor[0]
            
            cursor[1] = min_y-1
            
            pos_df = norm_pos_df

        # get the bounds of the positions
        min_x, min_y = pos_df.min(axis=1)
        max_x, max_y = pos_df.max(axis=1)       
        comp_bounds = [
            min(min_x,comp_bounds[0]),
            max(max_x,comp_bounds[1]),
            min(min_y,comp_bounds[2]),
            max(max_y,comp_bounds[3])
        ]

        comp_pos.append((subg,pos_df))

    node_size = 0.35
    
    comp_bounds[0] = -1
    comp_bounds[1] = max(10,comp_bounds[1])
    comp_bounds[2] -= node_size*2
    comp_bounds[3] += node_size*2
    comp_size = [
        comp_bounds[1]-comp_bounds[0],
        comp_bounds[3]-comp_bounds[2]
    ]
    fig_factor = np.array([.5])
    figsize = fig_factor*np.array(comp_size)

    # expand the figure right to see all the labels
    expand_x = 6
    figsize[0] *= expand_x
    
    f = None
    if ax is None :
        f, ax = subplots(figsize=fig_factor*figsize)
        ax.set_xlim(left=comp_bounds[0],right=comp_bounds[1]*expand_x)
        ax.axis('off')
    else :
        ax.set_xlim(left=comp_bounds[0],right=comp_bounds[1])
        
    ax.set_ylim(bottom=comp_bounds[2],top=comp_bounds[3])

    for subg, pos in comp_pos :
        nodes,_ = zip(*list(subg.nodes.items()))
        
        # draw edges
        for edge in subg.edges :
            n1 = pos[edge[0]]
            n2 = pos[edge[1]]
            ax.add_patch(PathPatch(Path([n1,n2])))
        for node in subg.nodes :
            ax.add_patch(Circle(pos[node],radius=node_size,facecolor=colormap.get(node,'orange')))
            ax.text(pos[node][0]+node_size+0.05, pos[node][1],node,verticalalignment='center',fontsize=12)

    return f

def jacc_thresh_scan(jacc) :

    n_comps = {}
    for i in arange(0.05,0.95,0.01) :
        jacc_thresh = jacc.copy()
        jacc_thresh[jacc_thresh<i] = 0
        if jacc_thresh.size == 0 :
            continue
        jacc_g = nx.from_numpy_matrix(jacc_thresh.values)
        comps = list(nx.connected_components(jacc_g))
        n_comps[i] = [
            len(comps),
            mean(list(map(len,comps))),
            max(list(map(len,comps)))
        ]
    df = pandas.DataFrame(n_comps,index=['num_components','mean_num_genes_per_comp','num_genes_largest_comp']).T
    if df.size != 0 :
        f, ax = subplots(figsize=(20,5))
        df.plot(ax=ax)
        #ax.axis('tight')
        ax.grid()
        ax.set_xticks(df.index)
        ax.set_xticklabels(['{:.2f}'.format(_) for _ in df.index],rotation=90)
    else :
        print('No significant results')
        
def cluster_gsea(gsea,jacc_thresh=0.25, colormap={}, ax=None) :

    jacc = compute_jaccard(gsea)

    #jacc_thresh_scan(jacc)

    jacc[jacc<jacc_thresh] = 0
    
    if jacc.sum().sum() == 0 :
        print('No surviving genesets to cluster')
        return None
    
    f = cluster_graph(jacc,colormap=colormap, ax=ax)

    return f

def union_leadingEdge(df) :
    edge_cols = [_ for _ in df.columns if 'leadingEdge' in _]
    leadingEdge = df[edge_cols].fillna('').sum(axis=1).apply(lambda x: ' '.join(set(x.split())))
    return leadingEdge