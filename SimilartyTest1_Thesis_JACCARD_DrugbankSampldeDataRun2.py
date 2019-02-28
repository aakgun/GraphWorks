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
conn = redis.Redis('localhost')

train_set = ["The sky is blue.", "The sun is bright."]  # Documents
test_set = ["The sun in the sky is bright."]  # Query
#stopWords = stopwords.words('english')

#vectorizer = CountVectorizer()
##print vectorizer
#transformer = TfidfTransformer()
##print transformer
#
#trainVectorizerArray = vectorizer.fit_transform(train_set).toarray()
#testVectorizerArray = vectorizer.transform(test_set).toarray()
#print ('Fit Vectorizer to train set', trainVectorizerArray)
#print ('Transform Vectorizer to test set', testVectorizerArray)
#
#transformer.fit(trainVectorizerArray)
#
#print (transformer.transform(trainVectorizerArray).toarray())
#
#transformer.fit(testVectorizerArray)
#
#tfidf = transformer.transform(testVectorizerArray)
#print (tfidf.todense())
#print ("\nSimilarity Score [*] ",cosine_similarity(trainVectorizerArray[0:2], trainVectorizerArray))

class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)

#authenticate("localhost:7474", "neo4j", "1234")
#graphneo = Graph()
#graphneo.delete_all()  


g = nx.Graph()
#g=nx.generators.directed.random_k_out_graph()
    
#object_drugbank,object,precipitant_drugbank,precipitant
with open("../PDDI-Datasets/DrugBank/drugbank4-DDIs-commaDelimetersample5.csv") as fcsv:
#with open("../PDDI-Datasets/DrugBank/drugbank4-DDIs-commaDelimeter_Correction5.csv") as fcsv:    
    
    i=1
    nodeid=0
    nodenr=0

    #dict = {'result':[{'key1':'value1','key2':'value2'}, {'key1':'value3','key2':'value4'}]}
    dict_nodes = {'nodes': []}
    dict_edges = {'links': []}
    
    dict_json= {"directed": False, "multigraph": False, "graph": {},'nodes': [],'links': [],'drugtitles': []} 
    dictnx_json= {'nxnodes': [],'nxlinks': []}
    dict_neighbours= {'neighbours': []} 
    dict_nodeneighbours= {'nodeneighbours': []} 
    
    OwlContent=""
    RDFContent=""
    RDFSContent=""
    RDFSContent2=""
    RDFSContentTemp=""
    OwlContent2=""
    OwlContent3=""
    OWLHeader=""
    
    RDFHeader = ("<?xml version=\"1.0\"?>\n" +
                 "<rdf:RDF \n" +
                 "xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n" +
                 "xmlns:interaction=\"http://www.ddi.com/interaction#\">\n")
    
    
    RDFSHeader = (""
                  "<?xml version=\"1.0\"?> \n"   
                  "<rdf:RDF \n" 
                  " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n" 
                  " xmlns:rdfs=\"http://www.w3.org/2000/01/rdf-schema#\" \n" 
                  " xmlns:foaf=\"http://xmlns.com/foaf/0.1/\" \n" +
                  " xmlns:owl=\"http://www.w3.org/2002/07/owl#\" \n" +
                  #" xmlns:owl=\"http://ddi.org/resource/interact/\" \n" +               
                  " xmlns:drugbank=\"http://bio2rdf.org/drugbank_vocabulary#\" \n" +                  
                  " xml:base=\"http://www.ddi.com/ddinteractions#\">\n")
                  
    RDFSHeader = (""                  
                  "@base <http://ddi.org/drugbank-resource/> . \n" +
                  "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> . \n" +
                  "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> . \n" +
                  "@prefix sesame: <http://www.openrdf.org/schema/sesame#> . \n" +
                  "@prefix owl: <http://www.w3.org/2002/07/owl#> . \n" +
                  "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .\n" +
                  "@prefix drugbank: \"http://bio2rdf.org/drugbank_vocabulary#\" \n" +                     
                  "@prefix fn: <http://www.w3.org/2005/xpath-functions#> . \n")
            
    OWLHeader=   """
<?xml version="1.0"?>
 <rdf:RDF
     xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
     xmlns:protege="http://protege.stanford.edu/plugins/owl/protege#"
     xmlns:obo="http://purl.obolibrary.org/obo#"
     xmlns:ace_lexicon="http://attempto.ifi.uzh.ch/ace_lexicon#"
     xmlns="http://purl.obolibrary.org/obo/DINTO_PKO_BRO.owl#"
     xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
     xmlns:purl="http://purl.obolibrary.org#"
     xmlns:obo2="http://purl.obolibrary.org/obo/"
     xmlns:owl="http://www.w3.org/2002/07/owl#"
     xmlns:skos="http://www.w3.org/2004/02/skos/core#"
     xmlns:iao="http://purl.obolibrary.org/obo/iao/"
     xmlns:biositemap="http://bioontology.org/ontologies/biositemap.owl#"
     xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
     xml:base="http://purl.obolibrary.org/obo/DINTO_PKO_BRO.owl" \n> 
                    """            


    
    OWLHeader3 =("<?xml version=\"1.0\"?>\n" +
                "<rdf:RDF xmlns=\"http://www.w3.org/2002/07/owl#\"\n" +
                "xml:base=\"http://www.w3.org/2002/07/owl\"\n" +
                "xmlns:void=\"http://rdfs.org/ns/void#\"\n" +
                "xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n" +
                "xmlns:terms=\"http://purl.org/dc/terms/\"\n" +
                "xmlns:owl=\"http://www.w3.org/2002/07/owl#\"\n" +
                "xmlns:xml=\"http://www.w3.org/XML/1998/namespace\"\n" +
                "xmlns:bio2rdf_vocabulary=\"http://bio2rdf.org/bio2rdf_vocabulary:\"\n" +
                "xmlns:xsd=\"http://www.w3.org/2001/XMLSchema#\"\n" +
                "xmlns:drugbank_vocabulary=\"http://bio2rdf.org/drugbank_vocabulary:\"\n" +
                "xmlns:rdfs=\"http://www.w3.org/2000/01/rdf-schema#\">\n" +
                "<Ontology/>\n"
                 "<ObjectProperty rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:ddi-interactor-in\"/>\n" +
                 "<ObjectProperty rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:ddi-interactor-out\"/>\n" +                 
                 "<Class rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Drug\"/>\n" +
                 "<Class rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Drug-Drug-Interaction\"/>\n" +
                 "<Class rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Resource\"/>\n" +
                 "<Class rdf:about=\"http://www.w3.org/2000/01/rdf-schema#Resource\"/>\n" +
                "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Drug\">\n" +
                "    <rdf:type rdf:resource=\"http://bio2rdf.org/drugbank_vocabulary:Resource\"/>\n" +
                "</NamedIndividual>\n" +
                "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Drug-Drug-Interaction\">\n" +
                "    <rdf:type rdf:resource=\"http://bio2rdf.org/drugbank_vocabulary:Resource\"/>\n" +
                "</NamedIndividual>\n" +
                "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:Resource\">\n" +
                "    <rdf:type rdf:resource=\"http://www.w3.org/2000/01/rdf-schema#Resource\"/>\n" +
                "</NamedIndividual>\n" +
                "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:ddi-interactor-in\">\n"
                "    <rdf:type rdf:resource=\"http://bio2rdf.org/drugbank_vocabulary:Resource\"/>\n" +
                "</NamedIndividual>  \n"
                "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank_vocabulary:ddi-interactor-out\">\n"
                "    <rdf:type rdf:resource=\"http://bio2rdf.org/drugbank_vocabulary:Resource\"/>\n" +
                "</NamedIndividual>  \n"
                 )
        
    i=0
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

        OwlContent3=OwlContent3 + ( 
                            "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank:"+drug2ID+"\">\n" +
                            "<drugbank_vocabulary:ddi-interactor-out rdf:resource=\"http://bio2rdf.org/drugbank_resource:"+drug1ID+"_"+drug2ID+"\"/>\n" +
                            "<drugbank_vocabulary:Drug rdf:resource=\"http://bio2rdf.org/drugbank/"+drug2ID+"\"/>\n" +                            
                            "</NamedIndividual>\n") 
        OwlContent3=OwlContent3 + ( 
                            "<NamedIndividual rdf:about=\"http://bio2rdf.org/drugbank:"+drug1ID+"\">\n" +
                            "<drugbank_vocabulary:ddi-interactor-in rdf:resource=\"http://bio2rdf.org/drugbank_resource:"+drug1ID+"_"+drug2ID+"\"/>\n" +
                            "<drugbank_vocabulary:Drug rdf:resource=\"http://bio2rdf.org/drugbank/"+drug1ID+"\"/>\n" +                            
                            "</NamedIndividual>\n")      
                
        append=True
        for item in dict_json["nodes"]:
                    #print(item)
                if item["id"] == drug1URI:
                        append=False
                        #print('x')
        if(append==True):
                    dict_json['nodes'].append(
                            {u"type": "class",u"id": drug1URI,u"nodeid": nodenr}
                            )
                    nodenr=nodenr+1 

        append=True
        for item in dict_json["nodes"]:
                if item["id"] == drug2URI:
                     append=False
                        #print('x')
        if(append==True):
                dict_json['nodes'].append(
                            {u"type": "class",u"id": drug2URI,u"nodeid": nodenr}
                            )
                nodenr=nodenr+1
                OwlContent=OwlContent + (                
                        "<NamedIndividual rdf:about=\"http://ddi.org/resource/drug/id/"+drug2ID +"\">\n" +
                        "<rdf:type rdf:resource=\"http://ddi.org/resource/Drug\"/>\n" +
                        #"<airport:tz>Asia/Manila</airport:tz>\n"  +
                        #"<rdfs:label>Ninoy Aquino Intl</rdfs:label>\n" +
                        "</NamedIndividual>\n")   
                OwlContent2=OwlContent2 + ( 
                        "<owl:Class rdf:about=\"http://purl.obolibrary.org/obo/"+drug2ID+"\"/>\n"
                        )
                    
       
        for item in dict_json['nodes']:
                    #print(item)
                    if item["id"] == drug1URI:
                        sourceid=item["nodeid"]
                        g.add_node(drug1URI, type='Class',weight=0,w2=0,w3=0,w4=0,w5=0,betweenness=0,degree_centrality=0,closeness_centrality=0,neighbors='')
                        #for nodex in g.nodes():
                        #    if (nodex == sourceid):
                        #            addnode=False
                        #if (addnode == True):
                        
        for item in dict_json['nodes']:
                    #print(item)
                    if item["id"] == drug2URI:
                        targetid=item["nodeid"]
                        g.add_node(drug2URI, type='Class',weight=0,w2=0,w3=0,w4=0,w5=0,betweenness=0,degree_centrality=0,closeness_centrality=0,neighbors='')
        ddi =   {u"type": "SubClassOf",
                         u"source": sourceid,
                         u"target": targetid,                 
                         u"drug1id": drug1URI,
                         u"drug2id": drug2URI,
                         u"object": obj, 
                         u"precipitant": pre
                         }
        
        dict_json["links"].append(ddi)
        g.add_edge(drug1URI,drug2URI, type='SubClassOf',weight=0,color='b')
fcsv.close

RDFEnd = ("</rdf:RDF>")     





with open("C:\public-PDDI-analysis-master\RDFXMLOWLDataSample.owl","w") as frdf:    
        #frdf.write(RDFHeader+RDFSContent+RDFEnd)
        frdf.write(OWLHeader3+OwlContent3+RDFEnd)
frdf.close()   

with open("../analysis-results/JsonNetworkxGraphFileOWLDataSample.json","w") as f:
    f.write(json.dumps(dict_json,cls=JSONEncoder, indent=4))
f.close()

with open("../analysis-results/JsonNetworkxGraphFileOWLDataSample.json") as fjson:
    json_data = json.load(fjson)
    fjson.close()   

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



dict_nodes = []    
print("Node ve Degreeleri Listeye Alma========================")
for n,nbrs in g.adj.items():
            n=n.strip()
#            print("Node:" + n  + " Degree: " + str(g.node[n]['weight'])) 
            dict_nodes.append(
                    {u"type": "class",u"nodeid": n,u"nodeweight": g.node[n]['weight'],u"nodeweight2": g.node[n]['w2']}
            )   



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

print("\n top 11 degree centrality connected component:")
#print ( nx.degree_centrality(g))
print("\n top 11 degree centrality connected component:")
for i in range(11):
    print(ordered_dcGC[i])
##print(dcGC_dict)


bb = nx.betweenness_centrality(g)
dc = nx.degree_centrality(g)
#nx.set_node_attributes(g, 'betweenness', bb)

##########################################################################
##print (bb)
sparql = SPARQLWrapper("http://localhost:3030/ds6/update")
#
for x in bb:
#    print(x +":" + str(bb.get(x)))
    g.node[x]['betweennes_centrality']=bb.get(x)
    sparqlquery = """
       PREFIX drugbank_vocabulary: <http://bio2rdf.org/drugbank_vocabulary:>
       delete {
          ?s ?p ?o 
       }
       where { 
          FILTER (?s = <http://bio2rdf.org/drugbank:""" + x.strip() + """> )
          FILTER (?p = drugbank_vocabulary:betweennes_centrality)
            ?s ?p ?o
       };        
       INSERT DATA
        { 
         <http://bio2rdf.org/drugbank:""" + x.strip() + """> drugbank_vocabulary:betweennes_centrality """ + str(bb.get(x)) + """  .
        }
    """
#    print(sparqlquery)
    sparql.setQuery(sparqlquery)
    sparql.query()

for x in dc:
#    print(x +":" + str(dc.get(x)))
    g.node[x]['degree_centrality']=dc.get(x)
    sparql.setQuery("""
       PREFIX drugbank_vocabulary: <http://bio2rdf.org/drugbank_vocabulary:>
       delete {
          ?s ?p ?o 
       }
       where { 
          FILTER (?s = <http://bio2rdf.org/drugbank:""" + x.strip() + """> )
          FILTER (?p = drugbank_vocabulary:degree_centrality)
            ?s ?p ?o
       };       
       INSERT DATA
        { 
         <http://bio2rdf.org/drugbank:""" + x.strip() + """> drugbank_vocabulary:degree_centrality """ + str(dc.get(x)) + """  .
        }
    """)
    sparql.query()

for n, nbrs in g.adj.items():

    sparqlquery = """
       PREFIX drugbank_vocabulary: <http://bio2rdf.org/drugbank_vocabulary:>
       delete {
          ?s ?p ?o 
       }
       where { 
          FILTER (?s = <http://bio2rdf.org/drugbank:""" + n.strip() + """> )
          FILTER (?p = drugbank_vocabulary:degree)
            ?s ?p ?o
       };     
       INSERT DATA
        { 
         <http://bio2rdf.org/drugbank:""" + n.strip() + """> drugbank_vocabulary:degree """ + str(g.degree(n)) + """  .
        }
    """
##    print(sparqlquery)
    sparql.setQuery(sparqlquery)
    sparql.query()
#########################################################################


print("Node Title Serviceden Alma ve Locale Set Etme")
sparql = SPARQLWrapper("http://localhost:3030/ds6/query") 
sparqlupdate = SPARQLWrapper("http://localhost:3030/ds6/update")   
data_corpus=[]
data_corpusid=[]

########################################################################
for n, nbrs in g.adj.items():
    n=n.strip()
 
    sparqlquery = """
            select 
                   distinct ?drug ?drugid ?drugdesc 
                    where {
                  ?drug <http://bio2rdf.org/drugbank_vocabulary:Drug> ?o
                    FILTER (?drug = <http://bio2rdf.org/drugbank:""" + n.strip() + """>)
                     SERVICE <http://bio2rdf.org/sparql>{
                       select *  where {
                  			?drug <http://purl.org/dc/terms/description> ?drugdesc.
                            ?drug <http://bio2rdf.org/bio2rdf_vocabulary:identifier> ?drugid
                  			FILTER (?drug = <http://bio2rdf.org/drugbank:""" + n.strip() +""">)
            			   }  
            	 		}
                  } 
        """
        
    sparqlquery = """
                select distinct ?drug ?drugdesc ?drugid    
                where {
                #?drug <http://bio2rdf.org/bio2rdf_vocabulary:identifier> '""" + n.strip() + """'
                			SERVICE <http://bio2rdf.org/sparql>{
                                       select *  where {
                                  			?drug <http://purl.org/dc/terms/description> ?drugdesc.
                                            ?drug <http://bio2rdf.org/bio2rdf_vocabulary:identifier> ?drugid
                                  			FILTER (?drug = <http://bio2rdf.org/drugbank:""" + n.strip() + """>)
                            			   }  
                            	 		}
                  } 
                     """  
    
   
#    print(sparqlquery)
    sparql.setQuery(sparqlquery)
#    sparql.query()        
    
    #
    # JSON example
    #print ('\n\n*** JSON Example')
    ret=""
    sparql.setReturnFormat(JSON)
    ret = sparql.query()
    results = ret.convert()    
    #print (results)
    #for result in sparql.query().bindings:
    #    print('%s' % (result["drugdesc"].value))
    
    drugdesc=""
    for result in results["results"]["bindings"]:
        drugid=result["drugid"]["value"]
        drugdesc=result["drugdesc"]["value"]
        drugdesc=drugdesc.replace("'","")
        g.node[n]['title']=drugdesc
        #print ("drug",result["drug"]["value"])
        #print ("drugdesc",result["drugdesc"]["value"]) 
        data_corpus.append(drugdesc)
        data_corpusid.append(drugid)
        dict_json['drugtitles'].append(
                            {u"drugid": result["drugid"]["value"],u"drugdesc": result["drugdesc"]["value"]}
                            )
   
    drugdesc=drugdesc.strip()
    drugdesc=drugdesc.replace('\r','\\r').replace('\n','\\n')
    sparqlqueryupdate = """
           PREFIX drugbank_vocabulary: <http://bio2rdf.org/drugbank_vocabulary:>
     
           delete {
              ?s ?p ?o 
           }
           where { 
              FILTER (?s = <http://bio2rdf.org/drugbank:""" + n.strip() + """> )
              FILTER (?p = drugbank_vocabulary:title)
                ?s ?p ?o
           };  
            
           INSERT DATA
            { 
             <http://bio2rdf.org/drugbank:""" + n.strip() + """> drugbank_vocabulary:title '""" + drugdesc + """'  .
            }
        """
#    print(sparqlqueryupdate)
    sparqlupdate.setQuery(sparqlqueryupdate)
#    sparqlupdate.query()

#########################################################################


#
#
sparql = SPARQLWrapper("http://localhost:3030/ds6/query")
###http://localhost:3030/ds6/query
sparql.setQuery("""
SELECT  ?drugid  ?drugdesc
WHERE {
 ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc
}
""")

data_corpus=[]
data_corpusid=[]

    #print (results)
sparql.setReturnFormat(JSON)
ret = sparql.query()
results = ret.convert()   
for result in results["results"]["bindings"]:
    drugid=result["drugid"]["value"]
    drugdesc=result["drugdesc"]["value"]
    drugdesc=drugdesc.replace("'","")
    g.node[drugid.split(":")[2]]['title']=drugdesc
    #print ("drug",result["drug"]["value"])
    #print ("drugdesc",result["drugdesc"]["value"]) 
    data_corpus.append(drugdesc)
    data_corpusid.append(drugid)
#    print(drugdesc)
    dict_json['drugtitles'].append(
                            {u"drugid": result["drugid"]["value"],u"drugdesc": result["drugdesc"]["value"]}
                            )

#########################################################################
#
#
#
#
#
#
### (?drug = <http://bio2rdf.org/drugbank:""" + n.strip() +""">)
##sparql.setReturnFormat(JSON)
##ret = sparql.query()
##results = ret.convert()   
##
##sparql.setReturnFormat(JSON)
##results = sparql.query().convert()
##
##
##print(results["results"]["bindings"])
##for result in results["results"]["bindings"]:
##    drugid=result["drugid"]["value"]
##    drugdesc=drugdesc + " " + result["drugdesc"]["value"]
##    print("drugdesc")    
##    
#
##pos=nx.spring_layout(g) # positions for all nodes
##nx.draw(g, pos, font_size=8, with_labels=True)
#
###import matplotlib.pyplot as plt
###plt.figure(figsize=(80, 60))
##pos=nx.spring_layout(g) # positions for all nodes
####edge_colors = [e[2]['color'] for e in g.edges(data=True)]
##nx.draw_networkx_edges(g,pos,edge_color='b',node_color='r',node_size=10   )
####nx.draw(g,pos,edge_color='b',node_color='r',node_size=10   )
##nx.draw(g, with_labels=True)
####nx.draw(g, with_labels=True)
####plt.savefig("path.png")                    
##nx.draw_networkx_nodes(g,pos,node_size=20)
#                    
#
##nx.draw(g, pos, font_size=8, with_labels=True)
##nx.spectral_layout(g)
##for p in pos:  # raise text positions
##    pos[p][1] += 0.1
##nx.draw_networkx_labels(g, pos,font_size=10,arrowstyle='->', arrowsize=10, width=2,arrows=True)
##nx.draw_networkx_edges(g, pos, arrowstyle='->', arrowsize=10, width=2)
##nx.convert_node_labels_to_integers(g)
##plt.savefig("../analysis-results/weighted_graph.png") # save as png
##plt.show()
#
##
dict_neighbors=[]
data_neighbors=[]
data_neighborsid=[]
for n, nbrs in g.adj.items():
    n=n.strip()
    g.node[n]['weight']=g.degree(n)
    dict_nodeneighbours['nodeneighbours'].append({u"nodeid":n,u"neighbours":list(g.neighbors(n))})
    tmp_neighbourwtotal=0
    tmp_neighbour2wtotal=0 
    dict_neighbors.append(
                    {u"node":  g.node[n],u"nodeid": n,u"nodeweight": g.node[n]['weight']}
            )    
    neighbourwtotal=0
    neighbour2wtotal=0    
    neighbors=""
    for x in g.neighbors(n):
        g.node[n]['w2']=g.degree(x)
        neighbors =  " " + x + neighbors
#        print("Node:" + n + " Degree: " + str(g.degree(n)) )
#        print("NodeNeighbour:" + x + " Degree: " + str(g.degree(x)) ) 
        neighbourwtotal = neighbourwtotal + g.degree(x)
        dict_neighbors.append(
                    {u"nodeneigh":  g.node[x],u"nodeneighid": x ,u"nodeneighw": g.node[n]['w2']}
            ) 
        for y in g.neighbors(x):
            neighbour2wtotal = neighbour2wtotal + g.degree(y)
#            print("NodeNeighbour2:" + y + " Degree: " + str(g.degree(y)) )
    g.node[n]['w2']=neighbourwtotal
    g.node[n]['w3']=neighbour2wtotal
    g.node[n]['neighbors']=neighbors
    data_neighbors.append(neighbors)
    data_neighborsid.append(n)
#
#
#tfidf_vectorizer = TfidfVectorizer()
#tfidf_matrix = tfidf_vectorizer.fit_transform(train_set)
#print (tfidf_matrix)
#cosine = cosine_similarity(tfidf_matrix[-1], tfidf_matrix)
#print (cosine)
#
#
#
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
ds = pd.read_csv("D://testSimalirty2.csv") #you can plug in your own list of products or movies or books here as csv file
tf = TfidfVectorizer(analyzer='word', ngram_range=(1, 3), min_df=0, stop_words='english')
#######ngram (1,3) can be explained as follows#####
##ngram(1,3) encompasses uni gram, bi gram and tri gram
##consider the sentence "The ball fell"
##ngram (1,3) would be the, ball, fell, the ball, ball fell, the ball fell
#
#
#import pandas as pd
#
df = pd.DataFrame(data_corpus,columns=['Title'])
df['ID'] = range(1, len(df) + 1)
df['NodeID'] = data_corpusid

#df = pd.DataFrame(ds,columns=['Title'])
#df['ID'] = range(1, len(df) + 1)

#
df2 = pd.DataFrame(data_neighbors,columns=['Neighbours'])
df2['ID'] = range(1, len(data_neighborsid) + 1)
df2['NodeID'] = data_neighborsid
#
#
#
tfidf_matrix = tf.fit_transform(ds['Title'])
tfidf_matrix2 = tf.fit_transform(df['Title'])
tfidf_matrix3 = tf.fit_transform(df2['Neighbours'])
cosine_similarities = cosine_similarity(tfidf_matrix,tfidf_matrix)
cosine_similarities2 = cosine_similarity(tfidf_matrix2,tfidf_matrix2)
cosine_similarities3 = cosine_similarity(tfidf_matrix3,tfidf_matrix3)




#
results = {} # dictionary created to store the result in a dictionary format (ID : (Score,item_id))
results2 = {} # dictionary created to store the result in a dictionary format (ID : (Score,item_id))
results3 = {} # dictionary created to store the result in a dictionary format (ID : (Score,item_id))
#
for idx, row in ds.iterrows(): #iterates through all the rows
#    # the below code 'similar_indice' stores similar ids based on cosine similarity. sorts them in ascending order. [:-5:-1] is then used so that the indices with most similarity are got. 0 means no similarity and 1 means perfect similarity
    similar_indices = cosine_similarities[idx].argsort()[:-5:-1] #stores 5 most similar books, you can change it as per your needs
    similar_items = [(cosine_similarities[idx][i], ds['ID'][i]) for i in similar_indices]
    results[row['ID']] = similar_items[1:]

#print("similar_indices",similar_indices)
#print("similar_items",similar_items)
#print("results",results)
#print("results",results)
#
#for idx, row  in df.iterrows(): #iterates through all the rows
##    print(idx)
##    print(row)
##    # the below code 'similar_indice' stores similar ids based on cosine similarity. sorts them in ascending order. [:-5:-1] is then used so that the indices with most similarity are got. 0 means no similarity and 1 means perfect similarity
##    similar_indices = cosine_similarities2[idx].argsort()[:-5:-1] #stores 5 most similar books, you can change it as per your needs
##    similar_items = [(cosine_similarities2[idx][i], df['Title'][i]) for i in similar_indices]
##    results2[row['Title']] = similar_items[1:]
#
#
for idx, row  in df.iterrows(): #iterates through all the rows
    # the below code 'similar_indice' stores similar ids based on cosine similarity. sorts them in ascending order. [:-5:-1] is then used so that the indices with most similarity are got. 0 means no similarity and 1 means perfect similarity
    similar_indices = cosine_similarities2[idx].argsort()[:-5:-1] #stores 5 most similar books, you can change it as per your needs
    similar_items = [(cosine_similarities2[idx][i], df['ID'][i]) for i in similar_indices]
    results2[row['ID']] = similar_items[1:]

for idx, row  in df2.iterrows(): #iterates through all the rows
    # the below code 'similar_indice' stores similar ids based on cosine similarity. sorts them in ascending order. [:-5:-1] is then used so that the indices with most similarity are got. 0 means no similarity and 1 means perfect similarity
    similar_indices = cosine_similarities3[idx].argsort()[:-6:-1] #stores 5 most similar books, you can change it as per your needs
    similar_items = [(cosine_similarities3[idx][i], df2['ID'][i]) for i in similar_indices]
    results3[row['ID']] = similar_items[1:]

#print(cosine_similarities[6].argsort()[:-5:-1])
#
#
##below code 'function item(id)' returns a row matching the id along with Book Title. Initially it is a dataframe, then we convert it to a list
def item(id):
#    return ds.loc[ds['ID'] == id]['Book Title'].tolist()[0]
    return df.loc[df['ID'] == id]['Title'].tolist()[0]
def itemNodeID(id):
#    return ds.loc[ds['ID'] == id]['Book Title'].tolist()[0]
    return df.loc[df['ID'] == id]['NodeID'].tolist()[0]

def recommend(id, num):
    if (num == 0):
        print("Unable to recommend any book as you have not chosen the number of book to be recommended")
    elif (num==1):
        print("Recommending " + str(num) + " book similar to " + item(id))
        
    else :
        print("Recommending " + str(num) + " books similar to " + item(id))
        
#    print("----------------------------------------------------------")
    recs = results2[id][:num]
    for rec in recs:
        print("You may also like to read: " + item(rec[1]) + " (score:" + str(rec[0]) + ")")


def similar(id):
    similarrec=[]        
    recs = results2[id][:5]
    for rec in recs:
        similarrec.append(itemNodeID(rec[1]))
#        print("Similars: " + item3NodeID(rec[1]) + " (score:" + str(rec[0]) + ")")
    return(similarrec)


#below code 'function item(id)' returns a row matching the id along with Book Title. Initially it is a dataframe, then we convert it to a list
def item3(id):
#    return ds.loc[ds['ID'] == id]['Book Title'].tolist()[0]
    return df2.loc[df2['ID'] == id]['Neighbours'].tolist()[0]
def item3NodeID(id):
#    return ds.loc[ds['ID'] == id]['Book Title'].tolist()[0]
    return df2.loc[df2['ID'] == id]['NodeID'].tolist()[0]
def recommend3(id, num):
    if (num == 0):
        print("Unable to recommend any book as you have not chosen the number of book to be recommended")
    elif (num==1):
        print("Recommending " + str(num) + " book similar to " + item3(id))
        
    else :
        print("Recommending " + str(num) + " books similar to " + item3(id))
        
    print("----------------------------------------------------------")
    recs = results3[id][:num]
    for rec in recs:
        print("You may also like to read: " + item3(rec[1]) + " (score:" + str(rec[0]) + ")")


def similar3(id):
    similarrec=[]            
    recs = results3[id][:5]
    for rec in recs:
        similarrec.append(item3NodeID(rec[1]))
#        print("Similars: " + item3NodeID(rec[1]) + " (score:" + str(rec[0]) + ")")
    return(similarrec)
        
#the first argument in the below function to be passed is the id of the book, second argument is the number of books you want to be recommended
print("----------------------------------------------------------")
recommend(3,5)    
print("----------------------------------------------------------")
recommend3(6,5)    

print(itemNodeID(3))
print(item3NodeID(6))

print(similar(3))
print(similar3(6))

#rval = json.dumps(results3)
#conn.set('results3', rval)
dfcosine = pd.DataFrame(cosine_similarities3)
conn.set('dfcosine', dfcosine.to_json(orient='split'))
#dfresults3 = pd.DataFrame(results3)
#conn.set('dfresults3', dfresults3.to_json(orient='split'))

#rval = json.dumps(json.dumps(results3, indent=4, sort_keys=True))
#conn.set('results3', rval)


print("-----------------Node Clusters-----------------------------------")
#print(nx.clustering(g))


#print("-----------------Short Paths-----------------------------------")
dict_nodesshortpath=[]

ErrorMsg=""
for eachnode in g.nodes():
    for eachnodex in g.nodes():
        if eachnodex != eachnode and g.degree(eachnode) > 0 and g.degree(eachnodex) > 0 :
            try:
#                print(eachnode,g.degree(eachnode))
#                print(eachnodex,g.degree(eachnodex))
                if len(nx.shortest_path(g,source=eachnode,target=eachnodex)) > 2 and  len(nx.shortest_path(g,source=eachnode,target=eachnodex)) < 4 :
                    dict_nodesshortpath.append({u"drugid1":eachnode,
                                            u"drugid2":eachnodex,
                                            u"pathlen":len(nx.shortest_path(g,source=eachnode,target=eachnodex)),
                                            u"pathlens":list(nx.shortest_path(g,source=eachnode,target=eachnodex))}
                )
            except nx.NetworkXNoPath:
                ErrorMsg="No path"
            
                
#            print("1: ", eachnode)
#            print("2: ", eachnodex)
#            print(nx.shortest_path(g,source=eachnode,target=eachnodex))    
            if len(list(nx.common_neighbors(g, eachnode, eachnodex))) > 0:
#                print("1: ", eachnode)
#                print("2: ", eachnodex)
                dict_neighbours['neighbours'].append({u"drugid1":eachnode,u"drugid2":eachnodex,u"neighbours":list(nx.common_neighbors(g, eachnode, eachnodex))})
#                dict_json['drugtitles'].append({u"drugid": result["drugid"]["value"],u"drugdesc": result["drugdesc"]["value"]})                
#            if len(list(nx.common_neighbors(g, eachnode, eachnodex))) > 0:
#                for x in nx.common_neighbors(g, eachnode, eachnodex):               
##                print("common nodes:" , nx.common_neighbors(g, eachnode, eachnodex))
#                    print(x)
#    print(g.node[eachnode].
print("----------------------------------------------------------")                
#print(dict_neighbours['neighbours'])          
#print("dict_nodesshortpath",dict_nodesshortpath)                
#
psearchDrug="DB00004"

#print(nx.strongly_connected_components(g))
dict_inpathsimilars=[]
dict_inpathsimilars2=[]
dict_inpathsimilars3=[]
orsimilarpathfilter=""
for eachdrugnode in dict_nodesshortpath:
#    print("eachdrugnode ",eachdrugnode,"--",eachdrugnode['pathlens'])
    for pathnodex in eachdrugnode['pathlens']:
#        print("pathnodex:",pathnodex)
        if psearchDrug == pathnodex:
#            print(pathnodex)
            dict_inpathsimilars.append({u"pathsimilarnode":pathnodex,u"pathsimilars":[eachdrugnode['drugid1'],eachdrugnode['drugid2']]}) 
            dict_inpathsimilars2.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid1']})
            dict_inpathsimilars2.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid2']})
            dict_inpathsimilars3.append(eachdrugnode['drugid1'])
            dict_inpathsimilars3.append(eachdrugnode['drugid2'])
#            dict_inpathsimilars.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid2']}) 
dict_inpathsimilars4 = set(dict_inpathsimilars3)
for i in dict_inpathsimilars4:
#    print (i)
    orsimilarpathfilter = orsimilarpathfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + i + ">" 
    orsimilarpathfilter = orsimilarpathfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + i + ">" 

print("---------------------orsimilarpathfilter-------------------------------------")        
#print(orsimilarpathfilter)
sparql = SPARQLWrapper("http://localhost:3030/ds6/query")


#print(list(df2['NodeID']))
ix=1
orsimilarfilter=""
for eachdrugnode in df2['NodeID']:
#    print (eachdrugnode, item3NodeID(ix))
    if (eachdrugnode == psearchDrug):
        for x in similar3(ix):
#            print(x)
            orsimilarfilter = orsimilarfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + x + ">" 
            orsimilarfilter = orsimilarfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + x + ">" 
        break
        
    ix=ix+1

print("---------------------orsimilarfilter-------------------------------------")            
#print(orsimilarfilter)
#    df2.loc[df2['ID'] == id]


a=1
orfilter=""
for eachnitem in dict_neighbours['neighbours']:
#    print(eachnitem['drugid1'],eachnitem['drugid2'],eachnitem['neighbours'],len(eachnitem['neighbours']))
    for eachnx in eachnitem['neighbours']:
        if (eachnx == psearchDrug):
            orfilter = orfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" +eachnitem['drugid1'] + ">" 
            orfilter = orfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" +eachnitem['drugid2'] + ">" 
            break
#            print (orfilter)
#        print(eachnx)
            a = a+1     
##    print(eachnitem.get("drugid1"))

tempquerystr = """                
SELECT distinct ?drugid  ?drugdesc
WHERE {
 ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc
 FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orfilter + orsimilarfilter + orsimilarpathfilter +""")
} limit 100
"""

tempquerystr2 = """       
SELECT * WHERE
{ 
{
    SELECT distinct ?drugid  ?drugdesc ("SearchItem" as ?Type)  ?betweennes ?degree  ?nodedegree
    WHERE {
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
      FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orfilter  +""")
    } order by desc( ?nodedegree)
}
UNION
{
     SELECT distinct ?drugid  ?drugdesc ("SimilarItem" as ?Type)  ?betweennes ?degree  ?nodedegree
     WHERE {
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
    FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orsimilarfilter  +""")
     } order by desc( ?nodedegree)
}
UNION
{
     SELECT distinct ?drugid  ?drugdesc ("PathSimilarItem" as ?Type)  ?betweennes ?degree  ?nodedegree
     WHERE {
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
       ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
     FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orsimilarpathfilter  +""")
     } order by desc( ?nodedegree)
}     
}
"""

print(tempquerystr)
print(tempquerystr2)
sparql.setQuery(tempquerystr)
sparql.setReturnFormat(JSON)
ret = sparql.query()
results = ret.convert()   
for result in results["results"]["bindings"]:
    print(result["drugid"]["value"],"----",result["drugdesc"]["value"],"\n")
    drugid=result["drugid"]["value"]
    drugdesc=result["drugdesc"]["value"]
    
def jaccard_similarity(query, document):
    intersection = set(query).intersection(set(document))
    union = set(query).union(set(document))
    return len(intersection)/len(union)

#print(data_neighbors)
data_neighbors.pop()
data_neighborsid.pop()
tokenize = lambda doc: doc.lower().split(" ")
all_documents = data_neighbors
tokenized_documents = [tokenize(d) for d in all_documents] # tokenized docs
all_tokens_set = set([item for sublist in tokenized_documents for item in sublist])

print(tokenized_documents[3])
print(tokenized_documents[6])
print(jaccard_similarity(tokenized_documents[3],tokenized_documents[6]))

#"from redisworks import Root
#root = Root()

rval = json.dumps(dict_nodesshortpath)
conn.set('key3', rval)

DemoQuery=tempquerystr2
print(DemoQuery)
    
def DrugSimilarSearchFunction (pDrugId):
    
    psearchDrug=pDrugId

#    root.something = dict_nodeneighbours
#    conn.hmset("pythonDict1", dict_nodeneighbours)

    rval = json.dumps(dict_nodesshortpath)
    conn.set('dict_nodesshortpath', rval)

    #print(nx.strongly_connected_components(g))
    dict_inpathsimilars=[]
    dict_inpathsimilars2=[]
    dict_inpathsimilars3=[]
    orsimilarpathfilter=""
    

        
    for eachdrugnode in dict_nodesshortpath:
    #    print("eachdrugnode ",eachdrugnode,"--",eachdrugnode['pathlens'])
        for pathnodex in eachdrugnode['pathlens']:
    #        print("pathnodex:",pathnodex)
            if psearchDrug == pathnodex:
    #            print(pathnodex)
                dict_inpathsimilars.append({u"pathsimilarnode":pathnodex,u"pathsimilars":[eachdrugnode['drugid1'],eachdrugnode['drugid2']]}) 
                dict_inpathsimilars2.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid1']})
                dict_inpathsimilars2.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid2']})
                dict_inpathsimilars3.append(eachdrugnode['drugid1'])
                dict_inpathsimilars3.append(eachdrugnode['drugid2'])
    #            dict_inpathsimilars.append({u"pathsimilarnode":pathnodex,u"pathsimilars":eachdrugnode['drugid2']}) 
    dict_inpathsimilars4 = set(dict_inpathsimilars3)
    for i in dict_inpathsimilars4:
    #    print (i)
        orsimilarpathfilter = orsimilarpathfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + i + ">" 
        orsimilarpathfilter = orsimilarpathfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + i + ">" 
    
    print("---------------------orsimilarpathfilter-------------------------------------")        
    #print(orsimilarpathfilter)
    sparql = SPARQLWrapper("http://localhost:3030/ds6/query")
    
#    df2list = df2.values.tolist()
#    rval = json.dumps(df2list)
#    conn.set('df2list', rval)


    conn.set('df2list', df2.to_json(orient='split'))

    #print(list(df2['NodeID']))
    
    
    ix=1
    orsimilarfilter=""
    
    for eachdrugnode in df2['NodeID']:
    #    print (eachdrugnode, item3NodeID(ix))
        if (eachdrugnode == psearchDrug):
            for x in similar3(ix):
    #            print(x)
                orsimilarfilter = orsimilarfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + x + ">" 
                orsimilarfilter = orsimilarfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" + x + ">" 
            break
            
        ix=ix+1
    
    print("---------------------orsimilarfilter-------------------------------------")            
    #print(orsimilarfilter)
    #    df2.loc[df2['ID'] == id]
    
    rval = json.dumps(dict_neighbours['neighbours'])
    conn.set('dict_neighbours', rval)
    
    a=1
    orfilter=""
    for eachnitem in dict_neighbours['neighbours']:
    #    print(eachnitem['drugid1'],eachnitem['drugid2'],eachnitem['neighbours'],len(eachnitem['neighbours']))
        for eachnx in eachnitem['neighbours']:
            if (eachnx == psearchDrug):
                orfilter = orfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" +eachnitem['drugid1'] + ">" 
                orfilter = orfilter + " " + " || ?drugid = <http://bio2rdf.org/drugbank:" +eachnitem['drugid2'] + ">" 
                break
    #            print (orfilter)
    #        print(eachnx)
                a = a+1     
    ##    print(eachnitem.get("drugid1"))
    
    tempquerystr = """                
    SELECT distinct ?drugid  ?drugdesc
    WHERE {
     ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc
     FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orfilter + orsimilarfilter + orsimilarpathfilter +""")
    } limit 100
    """
    
    tempquerystr2 = """       
    SELECT * WHERE
    { 
    {
        SELECT distinct ?drugid  ?drugdesc ("SearchItem" as ?Type)  ?betweennes ?degree  ?nodedegree
        WHERE {
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
          FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orfilter  +""")
        } order by desc( ?nodedegree)
    }
    UNION
    {
         SELECT distinct ?drugid  ?drugdesc ("SimilarItem" as ?Type)  ?betweennes ?degree  ?nodedegree
         WHERE {
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
        FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orsimilarfilter  +""")
         } order by desc( ?nodedegree)
    }
    UNION
    {
         SELECT distinct ?drugid  ?drugdesc ("PathSimilarItem" as ?Type)  ?betweennes ?degree  ?nodedegree
         WHERE {
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:title> ?drugdesc .
         ?drugid  <http://bio2rdf.org/drugbank_vocabulary:betweennes_centrality> ?betweennes .
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree_centrality> ?degree .
           ?drugid  <http://bio2rdf.org/drugbank_vocabulary:degree> ?nodedegree
         FILTER (?drugid = <http://bio2rdf.org/drugbank:""" + psearchDrug + """> """ + orsimilarpathfilter  +""")
         } order by desc( ?nodedegree)
    }     
    }
    """
    
    print(tempquerystr2)

    
