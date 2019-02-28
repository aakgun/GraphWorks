from __future__ import division
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import TfidfVectorizer
import numpy as np
import numpy.linalg as LA
from nltk.corpus import stopwords
import csv
import json
import networkx as nx
import string
import math
from py2neo import Node,Graph,authenticate
from SPARQLWrapper import SPARQLWrapper, JSON, XML, N3, RDF #, CSV, TSV
from SPARQLWrapper import SPARQLWrapper, XML, RDFXML, N3, JSONLD, JSON,  POST, GET, SELECT, CONSTRUCT
from networkx.readwrite import json_graph
import matplotlib.pyplot as plt
import redis
#conn = redis.Redis('localhost')



class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)

#authenticate("localhost:7474", "neo4j", "1234")
#graphneo = Graph()
#graphneo.delete_all()  


g = nx.DiGraph()
#g=nx.generators.directed.random_k_out_graph()
    
#object_drugbank,object,precipitant_drugbank,precipitant
with open("../../PDDI-Datasets/DrugBank/drugbank4-DDIs-commaDelimetersample5.csv") as fcsv:
#with open("../../PDDI-Datasets/DrugBank/drugbank4-DDIs-commaDelimeter_Correction5.csv") as fcsv:    
    
    nodeid=0
    nodenr=0

    #dict = {'result':[{'key1':'value1','key2':'value2'}, {'key1':'value3','key2':'value4'}]}
    dict_nodes = {'nodes': []}
    dict_edges = {'links': []}
    
    dict_json= {"directed": False, "multigraph": False, "graph": {},'nodes': [],'links': [],'drugtitles': []} 
    dictnx_json= {'nxnodes': [],'nxlinks': []}
    dict_neighbours= {'neighbours': []} 
    dict_nodeneighbours= {'nodeneighbours': []} 
    

    reader = csv.reader(fcsv, delimiter=',')
    for line in reader:  

        drug1ID=line[0].strip()
        drug2ID=line[2].strip()
        drug1URI="http://ddi.org/drugbank-resource/"+line[0].strip()
        drug2URI="http://ddi.org/drugbank-resource/"+line[2].strip()
        drug1URI=line[0].strip()
        drug2URI=line[2].strip()
        obj = line[1].strip()
        pre = line[3].strip()
        drug1id=int(line[0].strip()[-5:])
        drug2id=int(line[2].strip()[-5:])
        #print(drug1URI + '#'+ drug2URI)              
              
        drug =   {u"type": "class",
                         u"id": drug1URI
                         }      
        drug2 =  {u"type": "class",
                         u"id": drug2URI
                         }        
        ddi =   {u"type": "SubClassOf",
                         u"source": drug1URI,
                         u"target": drug2URI,                 
                         u"drug1id": drug1URI,
                         u"drug2id": drug2URI,
                         u"object": obj, 
                         u"precipitant": pre,  
                         }
  
#        g.add_node(drug1id)
#        g.add_edge(drug1id,drug2id)        
#        g.add_node(drug1id, type='Class',weight=0,w2=0,w3=0,w4=0,w5=0,betweenness=0,degree_centrality=0,closeness_centrality=0,neighbors='')
#        g.add_node(drug1URI, type='Class',weight=0,w2=0,w3=0,w4=0,w5=0,betweenness=0,degree_centrality=0,closeness_centrality=0,neighbors='')
#        g.add_edge(drug1id,drug2id, type='SubClassOf',weight=0,color='b')
        g.add_node(drug1URI,weight=0)
        g.add_edge(drug1URI,drug2URI)        
#       
#        g.add_edge(drug1URI,drug2URI, type='SubClassOf',weight=0,color='b')
fcsv.close

A = nx.to_numpy_matrix(g)

print (A)


import networkx as nx
from node2vec import Node2Vec
node2vec = Node2Vec(g, dimensions=64, walk_length=30, num_walks=200, workers=4)  # Use temp_folder for big graphs


# Precompute probabilities and generate walks - **ON WINDOWS ONLY WORKS WITH workers=1**
node2vec = Node2Vec(g, dimensions=64, walk_length=30, num_walks=200, workers=4)  # Use temp_folder for big graphs
# Embed nodes
model = node2vec.fit(window=10, min_count=1, batch_words=4)  # Any keywords acceptable by gensim.Word2Vec can be passed, `diemnsions` and `workers` are automatically passed (from the Node2Vec constructor)
# Look for most similar nodes
#model.wv.most_similar('DB00005')  # Output node names are always strings
## Save embeddings for later use
#model.wv.save_word2vec_format('EMBEDDING_FILENAME')
## Save model for later use
#model.save('EMBEDDING_MODEL_FILENAME')
## Embed edges using Hadamard method
#from node2vec.edges import HadamardEmbedder
#edges_embs = HadamardEmbedder(keyed_vectors=model.wv)
## Look for embeddings on the fly - here we pass normal tuples
#edges_embs[('1', '2')]
#''' OUTPUT
#array([ 5.75068220e-03, -1.10937878e-02,  3.76693785e-01,  2.69105062e-02,
#       ... ... ....
#       ..................................................................],
#      dtype=float32)
#'''
## Get all edges in a separate KeyedVectors instance - use with caution could be huge for big networks
#edges_kv = edges_embs.as_keyed_vectors()
## Look for most similar edges - this time tuples must be sorted and as str
#edges_kv.most_similar(str(('1', '2')))
## Save embeddings for later use
#edges_kv.save_word2vec_format('EDGES_EMBEDDING_FILENAME')

#######################################Neo4J#########################
query = """
        MATCH (n)
        DETACH DELETE n
        """
#graphneo.cypher.execute(query,json_data=json_data)  
print("query:",query)
query = """
       WITH {json_data} AS document
       UNWIND document.nodes AS node
       MERGE (:NodeId {name: node.id})
       """
#print("query2:",query)
#graphneo.cypher.execute(query,json_data=json_data)
i=0
for item in dict_json['links']:
        query = ("WITH {json_data} AS document \n" +
             " MATCH (a:NodeId),(b:NodeId) \n"  +
             " WHERE a.name=\'" +  str(item["drug1id"]) + "\'  AND b.name= \'" + str(item["drug2id"]) + "\' \n" + 
             " MERGE (a)-[r:DDI]->(b) \n" )
        i=i+1         
        #print (str(item["drug1id"]))
        #print("query3:",i,query)
#        graphneo.cypher.execute(query,json_data=json_data) 
#######################################Neo4J#########################

print("Degree Centrality")
#print(nx.degree_centrality(g))
#print(sorted(dict_nodes, key=itemgetter('nodeweight'),reverse=True))

dcGC_dict = nx.degree_centrality(g)
ordered_dcGC = sorted(dcGC_dict, key = dcGC_dict.get,reverse = True)

print('# of edges: {}'.format(g.number_of_edges()))
print('# of nodes: {}'.format(g.number_of_nodes()))

      
#indcGC_dict = nx.in_degree_centrality(g)
#ordered_indcGC = sorted(indcGC_dict, key = indcGC_dict.get,reverse = True)
#
#outdcGC_dict = nx.out_degree_centrality(g)
#ordered_outdcGC = sorted(outdcGC_dict, key = outdcGC_dict.get,reverse = True)

#print("\n top 11 degree centrality connected component:")
##print ( nx.degree_centrality(g))
#print("\n top 11 degree centrality connected component:")
#for i in range(11):
#    print(ordered_dcGC[i])
###print(dcGC_dict)


bb = nx.betweenness_centrality(g)
dc = nx.degree_centrality(g)
#nx.set_node_attributes(g, 'betweenness', bb)

#import matplotlib.pyplot as plt
#plt.figure(figsize=(80, 60))
pos=nx.spring_layout(g) # positions for all nodes
##edge_colors = [e[2]['color'] for e in g.edges(data=True)]
nx.draw_networkx_edges(g,pos,edge_color='b',node_color='r',node_size=10   )
##nx.draw(g,pos,edge_color='b',node_color='r',node_size=10   )
nx.draw(g, with_labels=True)
##nx.draw(g, with_labels=True)
##plt.savefig("path.png")                    
nx.draw_networkx_nodes(g,pos,node_size=20)


#plt.subplot(122)
#nx.draw(g, with_labels=True, font_weight='bold')
