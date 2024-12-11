Introducción al Análisis de redes
================
Guilermo de Anda-Jáuregui

**Hacking 4 Biopunks**

# Análisis de redes

# ¿Qué es una red?

Una red es un objeto matemático que tiene dos conjuntos:

  - Un conjunto de elementos, representados por *nodos*
  - Un conjunto de relaciones entre los elementos, representados por
    *enlaces*

# ![](redes_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# ¿Qué haremos hoy?

1 Leer una red (y escribirla, de una vez)

2 Describir la red: tamaño, caminos mas cortos, número de componentes,
etc.

3 Estudiar propiedades de los nodos: centralidades, coeficientes de
agrupamiento locales, etc.

4 Estudiar propiedades de los enlaces: intermediacion

5 Buscar grupos (comunidades)

6 Graficar



# Iniciemos

carguemos los paquetes


**R** 
``` r
library(igraph)
```

**Python** 
``` python
import networkx as nx
```

# Leer y escribir redes

Para el ejemplo, voy a usar la red de interacciones de los personajes en la novela de *Los
Miserables* (crédito: Donald Knuth, Mark Newman). 

Pueden seguir los ejercicios con cualquiera de las redes biológicas que les compartimos.

(tambien está la red de les Mis, si quieren replicar este código)

*Hay que descargar la red guardando el archivo haciendo click en "raw" y guardando la página. Nota: en algunos sistemas operativos / navegadores hay problema para descargar el raw file. Si no pueden descargarlo, intenten copiar y pegar el texto en un editor de textos y guardar el archivo.*



**Colocar el archivo en el directorio de trabajo**

## Leamos la red usando read.graph

**R** 
``` r
g <- read.graph(file = "les_miserables.graphml", format = "graphml")
```

**Python** 
``` python
g = nx.read_graphml("les_miserables.graphml")

```

Considerar que existen varios tipos de formatos de redes.

Examinemos nuestra red

**R**
``` r
g
```

    ## IGRAPH e2a34c4 U--- 77 254 -- 
    ## + attr: id (v/n), label (v/c)
    ## + edges from e2a34c4:
    ##  [1]  1-- 2  1-- 3  1-- 4  3-- 4  1-- 5  1-- 6  1-- 7  1-- 8  1-- 9  1--10
    ## [11] 11--12  4--12  3--12  1--12 12--13 12--14 12--15 12--16 17--18 17--19
    ## [21] 18--19 17--20 18--20 19--20 17--21 18--21 19--21 20--21 17--22 18--22
    ## [31] 19--22 20--22 21--22 17--23 18--23 19--23 20--23 21--23 22--23 17--24
    ## [41] 18--24 19--24 20--24 21--24 22--24 23--24 13--24 12--24 24--25 12--25
    ## [51] 25--26 24--26 12--26 25--27 12--27 17--27 26--27 12--28 24--28 26--28
    ## [61] 25--28 27--28 12--29 28--29 24--30 28--30 12--30 24--31 31--32 12--32
    ## [71] 24--32 28--32 12--33 12--34 28--34 12--35 30--35 12--36 35--36 30--36
    ## + ... omitted several edges

**Python** 
``` python
g
len(g.nodes()) # num. nodos
len(g.edges()) # num. aristas
```

Nuestra red tiene 77 nodos, con 254 enlaces entre ellos. Es una red *no
dirigida* y *no pesada*.


¿Cómo es *tu* red?


## Hagamos una primera visualización

``` r
plot(g)
```

**Python** 
``` python
nx.draw(g)
```

![](redes_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Esta visualización no es muy estética… ya la mejoraremos

## Exportemos la red

Existen varios formatos para representar redes: matriz de adyacencia, lista de enlaces, etc. 

Escribamos nuestra red en graphml, un formato basado en xml que permite portar información de nodos y aristas.

**R**
``` r
#escribir
write.graph(graph = g, 
            file = "red_mis.graphml", 
            format = "graphml"
            )
```

**Python**
``` python
#escribir
nx.write_graphml(file="red_mis.graphml")
```

### Para explorar: 

- Revisen el archivo de salida; ¿cómo se ve? 
- Experimenten otro tipo de archivos de salida: gml, edge list, json... ¿cómo difieren? 
- Para los usuarios de R y Python: ¿qué pasa si intento leer una red escrita en Python en R o viceversa? 


### Interludio: ¿Hay una red en mis datos? 

Respuesta: **es más probable de lo que te imaginas** 


Ejercicio opcional: 

- Descargar estos datos:

<https://raw.githubusercontent.com/guillermodeandajauregui/useR2019INMEGEN/master/movies.csv>


Se ven más o menos así: 

``` r
    ##               person      movie billing  role
    ## 1    Viggo Mortensen green book       1 actor
    ## 2     Mahershala Ali green book       2 actor
    ## 3   Linda Cardellini green book       3 actor
    ## 4 Dimiter D. Marinov green book       4 actor
    ## 5        Mike Hatton green book       5 actor
    ## 6        Iqbal Theba green book       6 actor
``` 

Preguntas:

- ¿Puedo leer esto como una red? (Pista: con Python, *read_edgelist()*, con R, leer la tabla y *graph_from_dataframe()*)
- ¿Qué tipo de red sería esta? 
- ¿Tienes algún conjunto de datos que se parezca a esta tabla?


Dejemos esta red para otro momento…

Recordar: si tenemos un data set con dos columnas que puedan
representarse como factores, en principio tenemos alguna estructura de
relaciones que podemos representar como una red. ¿Es una red
interesante? *No lo sabemos*.

### Fin del interludio 

Volvamos a nuestra red biológica (o literaria)...

## Analicemos propiedades globales de la red

  - Componentes

**R**
``` r
mis_componentes <- components(g)
```

**Python**

``` python
nx.number_connected_components(g)
nx.connected_components(g)
#para redes dirigidas hay que usar weak or strong connected components
nx.strongly_connected_components(g)
nx.weakly_connected_components(g)
```


- Longitud promedio de caminos más cortos

**R**
``` r
mis_caminos <- average.path.length(g)
```

**Python**
``` python
mis_caminos = nx.average_shortest_path_length(g)
``` 

- clustering coefficient global

**R**
``` r
mi_clustering <- transitivity(g, type = "global")
```

**Python**
``` python
average_clustering(g)
``` 

- Diámetro 

**R**
``` r
mi_diametro <- diameter(g)
mi_diametro
```

**Python**
``` python
nx.diameter(g)
``` 


## analizar propiedades de los nodos

Accedamos a los nodos


**R**
``` r
V(g)
```

**python**
``` python
g.nodes()
g.nodes(data=True)
```

Podemos tener la información completa de los nodos

**R**
``` r
mi_df_nodos <- get.data.frame(x = g, 
                              what = "vertices")
```

**Python**
``` python
#usando pandas
import pandas as pd
pd.DataFrame.from_dict(dict(g.nodes(data=True)), orient='index')
```


Podemos asignar nuevas propiedades a los nodos

**R**
``` r
#hagamos una nueva propiedad, name, que sea una copia de label
V(g)$name <- V(g)$label

mi_df_nodos <- get.data.frame(x = g, 
                              what = "vertices")

```

**Python**
``` python
my_labels = nx.get_node_attributes(g, "label")
nx.set_node_attributes(g, values=my_labels, name="name")
pd.DataFrame.from_dict(dict(g.nodes(data=True)), orient='index')
```

### Calculemos algunas medidas de centralidad

  - grado

**R**
``` r
V(g)$grado <- degree(g)
```

**Python**
``` python
nx.degeee(g)

```

  - intermediación (*betweenness*)


**R**
``` r
V(g)$intermediacion <- betweenness(g)
```

**Python**
``` python
nx.betweenness_centrality(g)
```

### Podemos tener también el coeficiente de agrupamiento local

**R**
``` r
V(g)$agrupamiento <- transitivity(g, type = "local", isolates = "zero")
```

**Python**
``` r
nx.clustering(g)
```

### Ejercicio opcional: 

¿Con lo que hemos aprendido, cómo obtendríamos la distribución de grado de nuestra red? ¿Qué nos dice sobre nuestra red?



## analizar propiedades de los enlaces

Accedamos a los enlaces

**R**
``` r
E(g)
```

**Python**
``` r
g.edges()
```

Definamos la propiedad de edge betweenness 

**R**
``` r
E(g)$intermediacion <- edge.betweenness(g)
```

**Python**
``` r
nx.edge_betweenness(g)
```

# Y ahora - un super poder para ustedes

## Clustering –\> Comunidades

Podemos usar la estructura de las redes para buscar un tipo especial de
clusters llamados *comunidades*. Estas son, intuitivamente, conjuntos de
nodos que comparten más enlaces entre ellos, que con nodos fuera del
conjunto.

Podemos usar diferentes algoritmos para detectar comunidades.

``` r
comm.louvain <- cluster_louvain(g)
comm.louvain
```

``` r
comm.louvain <- cluster_louvain(g)
comm.louvain
```

``` r

#podemos agregar la comunidad a la que pertenecen a cada nodo 
V(g)$comm.louvain <- membership(comm.louvain)
```

Tenemos varios algoritmos implementados

``` r
comm.infomap <- cluster_infomap(g)
comm.walktrp <- cluster_walktrap(g)
comm.fstgred <- cluster_fast_greedy(g)

V(g)$comm.infomap <- membership(comm.infomap)
V(g)$comm.walktrp <- membership(comm.walktrp)
V(g)$comm.fstgred <- membership(comm.fstgred)
```

Extraigamos una vez más la tabla de nodos

``` r
mi_df_nodos <- get.data.frame(x = g, 
                              what = "vertices")
mi_df_nodos %>% head()
```

## Hagamos unos plots

``` r
plot(g, 
     vertex.label.cex = 0.75,
     vertex.size = degree(g), 
     layout = layout_nicely
     )
```

Genera una visualización con el tamaño de los nodos proporcional al
grado, escogiendo heurísticamente el “mejor acomodo”, y con el tamaño de
las etiquetas más chico

``` r
plot(g, 
     vertex.label.cex = 0.75, 
     edge.curved = TRUE, 
     layout = layout_in_circle
     )
```

## Esta visualización acomoda los nodos en un círculo, y pone los enlaces curveados


### Protip para la gente de R: make it Tidy!

``` r
library(tidygraph)
library(ggraph)

gt <- as_tbl_graph(g)

gt %>% 
  ggraph(layout = 'kk') + 
  geom_edge_link(aes(), show.legend = FALSE) + 
  geom_node_point(aes(colour = as.factor(comm.louvain)), size = 7) + 
  guides(colour=guide_legend(title="Louvain community")) +
  theme_graph()
```
