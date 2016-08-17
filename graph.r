filter.large.cluster<-function(graph,min.n=3){
	commu<-walktrap.community(graph)
	which(crossing(commu,graph)==TRUE)->ind
	delete.edges(graph,ind)->tmp.graph
	decompose.graph(tmp.graph)->tmp.cluster.graph
	which(lapply(tmp.cluster.graph,vcount)>=min.n)->ind
	tmp.cluster.graph[ind]->tmp.cluster.graph
	return (tmp.cluster.graph)
}



filter.sparse.cluster<-function(cluster.graph,meth.info){
new.cluster<-NULL
for (i in 1:length(cluster.graph)){
cpg.name<-V(cluster.graph[[i]])$name
cpg.pos<-meth.info[cpg.name,]
cpg.cluster<-clusterMaker(cpg.pos$chr,cpg.pos$position,maxGap=2000)
cpg.pos$index<-cpg.cluster
if(sum(table(cpg.cluster) >= 3)>=1){
	append(new.cluster,cluster.graph[i])->new.cluster
}
}
return(new.cluster)
}

decompose.large.cluster<-function(cluster.graph,min.graph.density=0.3){
new.cluster.graph<-NULL
for (i in 1:length(cluster.graph)){
if (graph.density(cluster.graph[[i]])<min.graph.density){
	filter.large.cluster(cluster.graph[[i]])->filter.cluster
	append(new.cluster.graph,filter.cluster)->new.cluster.graph
}
else{
	#print (i)
	append(new.cluster.graph,cluster.graph[i])->new.cluster.graph
}
}

return (new.cluster.graph)
}
