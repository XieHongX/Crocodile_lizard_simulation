# Script to make SLIM job script
# USAGE: ./make_slim_bott.sh [Na] [Nb] [nF]
# This script is modified from Kyriazis's simulation work on Evolution Letters, 2020

# Set Na, the ancestral population size for burn-in
Na=${1}

# Set Nb, the final bottleneck population size
Nb=${2}

# Set nF, the number of founders:same as Nb in this case
nF=${3}

# Make script
cat > slim_bottleneck_gradual_contraction_age1_${Na}Na_${Nb}Nb_${nF}nF.slim << EOM

initialize() {
	
	initializeSLiMModelType("nonWF");
	defineConstant("K1", ${Na});
	defineConstant("K3", ${Nb});
	defineConstant("num_founders", ${nF});
	defineConstant("sampleSize", 30);
	defineConstant("g",20000); //number of genes
	defineConstant("ROHcutoff", 1000000);
	defineConstant("geneLength", 1500);
	defineConstant("seqLength", g*geneLength);
	//cat("Genome length:"+seqLength+"\n");	
	
        // half of all deleterious mutations (m1 or m2) become neutral here
        // so need to increase mutation rate by (2.31/3.31)/0.5=1.396 
        // to recover original volume of deleterious mutations 
        // this approach leads to more neutral mutations than the original model
        // so we dont draw those separately
        initializeMutationRate(4.258e-9);
        defineConstant("h_strDel", 0);
        defineConstant("h_wkModDel", 0.25);
        initializeMutationType("m1", h_strDel, "s", "x=rgamma(1,-0.01314833,0.186); if(x < -0.01){;return(x);}else{;return(0);}");
        initializeMutationType("m2", h_wkModDel, "s", "x=rgamma(1,-0.01314833,0.186); if(x >= -0.01){;return(x);}else{;return(0);}");
        initializeGenomicElementType("g1", c(m1,m2), c(1,1));	
	
	// approach for setting up genes on different chromosomes adopted from Jacqueline's wolf scripts 
	
	// vector of # genes on 16 different crocodile lizard chromosomes - need to rescale according to number of desired genes
	gene_nums=c(325,313,300,300,229,146,96,31,41,41,38,36,36,29,28,11);
	gene_nums = gene_nums*g/2000; //need to scale to number of desired genes since above array was originally set up for 2000 genes
	
	
	for (i in 1:g){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
		
	
	rates=NULL;
	
	// Multiple chromosomes:
	for (i in 1:(size(gene_nums)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[size(gene_nums)-1]-1)));
	
	ends=NULL;
	for (i in 1:g){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	
	initializeRecombinationRate(rates, ends);

}



// define function getStats that randomly samples a subpopulation for sampSize # of inds and outputs a string of: 
// pop size, mean fitness, heterozygosity, mean Froh, and avg num of variants of different classes per individual (very str del, str del, mod del, wk del)

function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F);
	
	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations??
	
	//initialize vectors
	ROH_length_sumPerInd = c();
	Num_VstrDel_muts = c();
	Num_strDel_muts = c();
	Num_modDel_muts = c();
	Num_wkDel_muts = c();
	ind_het = c();
	fitness_population = c();
	
	for (individual in i) {
		
		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed 
		indm_uniq = unique(indm);
		
		genotype = sapply(indm_uniq, "sum(indm == applyValue);");
		
		// tally number of mutations for different classes of selection coefficient per individual
		s = individual.genomes.mutations.selectionCoeff;
		
		Num_VstrDel_muts = c(Num_VstrDel_muts, sum(s<=-0.05));
		Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
		Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
		Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));
		
		if (isNULL(genotype)) {
			ind_het = c(ind_het, 0); //putting this here to avoid error when trying to sum null vector
			next;
		}
		
		ind_het = c(ind_het, sum(genotype==1)/(seqLength));
		
		//code for getting ROHs
		
		ID_het = (genotype == 1); //outputs T/F for genotypes if they are het or homDer
		ID_homDer = (genotype == 2);
		pos_het = indm_uniq.position[ID_het]; //outputs positions of heterozgoys genotypes
			
		startpos = c(0, pos_het); //adds 0 to beggining of vector of hets
		endpos = c(pos_het, sim.chromosome.lastPosition); //adds last position in genome to vector of hets
		pos_het_diff = endpos - startpos;
		ROH_startpos = startpos[pos_het_diff > ROHcutoff]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos = endpos[pos_het_diff > ROHcutoff];
		ROH_length = pos_het_diff[pos_het_diff > ROHcutoff]; //vector of ROHs for each individual	
		ROH_length_sum = sum(ROH_length);
		ROH_length_sumPerInd = c(ROH_length_sumPerInd, ROH_length_sum); // add sum of ROHs for each individual to vector of ROHs for all individuals
		
		// calculate individual fitness - code from Bernard	
		allmuts = c(individual.genomes[0].mutationsOfType(m1), individual.genomes[1].mutationsOfType(m1));
		uniquemuts = individual.uniqueMutationsOfType(m1);
		
                fitness_individual = c();
                
                if (size(uniquemuts) > 0){
                        for (u in uniquemuts){
                                places = (allmuts.id == u.id);
                                uu = allmuts[places];
                                if (size(uu) == 2) {
                                        fitness = 1 + sum(uu.selectionCoeff)/2;
                                } else if (size(uu) == 1) {
                                        if (u.mutationType == m1){
                                                fitness = 1 + uu.selectionCoeff * h_strDel;
                                        }
                                        if (u.mutationType == m2){
                                                fitness = 1 + uu.selectionCoeff * h_wkModDel;
                                        }
                                }
                                fitness_individual = c(fitness_individual, fitness);
                        }
                        fitness_individual = product(fitness_individual);
                        fitness_population = c(fitness_population, fitness_individual);
                } else {
                        fitness_population = c(fitness_population, 1);
                }
        }	
	return(pop.individuals.size() + "," + mean(fitness_population) + "," + mean(ind_het) + "," + mean(ROH_length_sumPerInd)/seqLength + "," + mean(Num_VstrDel_muts) + "," + mean(Num_strDel_muts)+ "," + mean(Num_modDel_muts) + "," + mean(Num_wkDel_muts));
}



reproduction() {
	if (individual.age > 1){
                mate = subpop.sampleIndividuals(1, minAge=2);   
                subpop.addCrossed(individual, mate);
        }
}



1 early() {
	cat("gen,popSize,meanFitness,meanHet,FROH,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", 10);
}

1:$((${Na}*10)) early() {
	p1.fitnessScaling = K1 / p1.individualCount;
}

//track statistics pre-bottleneck every 1000 generations
1:$((${Na}*10)) late() {
		if (sim.generation % 1000 == 0) {
			stats = getStats(p1, sampleSize);
			cat(sim.generation + "," + stats + "\n");
		}
}



// bottleneck to K=1000 for 1000 gens
$((${Na}*10+1)) early(){

        sim.addSubpop("p2",0);
        migrants = sample(p1.individuals, 1000);
        p2.takeMigrants(migrants);

}


// fitness scaling for p2

$((${Na}*10+1)):$((${Na}*10+1000)) early() {
	p1.fitnessScaling = 0; // kill off p1
	p2.fitnessScaling = 1000 / p2.individualCount;

}


//track statistics for intermediate contraction every 50 generations
$((${Na}*10+1)):$((${Na}*10+1000)) late() {
                if (sim.generation % 50 == 0) {
                        stats = getStats(p2, sampleSize);
                        cat(sim.generation + "," + stats + "\n");
                }
}




// bottleneck to p3
$((${Na}*10+1001)) early(){
	sim.addSubpop("p3",0);
	migrants = sample(p2.individuals, num_founders);
	p3.takeMigrants(migrants);
	cat("gen,K3,p_death,popSizeP3,meanFitness,meanHet,FROH,avgVStrDel,avgStrDel,avgModDel,avgWkDel,fixedVStrDel,fixedStrDel,fixedModDel,fixedWkDel" + "\n");
	
	sim.tag = K3; // use sim.tag to keep track of K3 from one generation to the next
}



// fitness scaling for p3  

$((${Na}*10+1001)):$((${Na}*10+6000)) early() {
	p2.fitnessScaling = 0; // kill off p2
	
	// kill off individuals at random - not sure if I should then adjust the individualCount
	inds = p3.individuals;
	
	//simulate beta distribution
	alpha = 0.5;
	beta = 8;
	x1 = rgamma(1, mean = alpha, shape=alpha);
	x2 = rgamma(1, mean = beta, shape=beta);
	beta = x1/(x1+x2); //probability of stochastic mortality this generation
	
	
	//set probability of death for each generation equal to outcome of beta 	
	for(i in inds){
		kill = rbinom(1,1,beta);
		if(kill==1){
			i.fitnessScaling = 0.0;
		}
	}
	
	//without environmental stochasticity
	sim.tag = K3;
	p3.fitnessScaling = sim.tag / p3.individualCount;
	
	cat(sim.generation + "," + sim.tag + "," + beta + ",");

}



// track statistics for P3 every generation and terminate when the population goes to 1 individual or after 5000 generations
$((${Na}*10+1001)):$((${Na}*10+6000)) late() {
	if(p3.individuals.size() < 2){
		stats_P3 = c("NA,NA,NA,NA,NA,NA,NA,NA"); //cant get stats from just one individual
	}
	if(p3.individuals.size() < sampleSize & p3.individuals.size() > 1){	// case when p3 size is less than sample size but greater than 1
		stats_P3 = getStats(p3, p3.individuals.size());
	}
	if(p3.individuals.size() >= sampleSize){ //case when p3 size is greater than or equal to sample size
		stats_P3 = getStats(p3, sampleSize);
	}
		
	//calculate fixed mutations in different classes
	 muts = sim.mutations;

	counts = sim.mutationCounts(NULL, muts);//a vctor of mutation counts for each mutation

	max = sum(sim.subpopulations.individualCount) * 2;
	s = sim.mutations.selectionCoeff;//a vector of selection coefficients
	fixed = s[counts == max];//all the fixed mutations
	fixed_VstrDel_muts = sum(fixed <=-0.05);
	fixed_strDel_muts = sum(fixed <=-0.01);
	fixed_modDel_muts = sum(fixed <=-0.001 & fixed > -0.01);
	fixed_wkDel_muts = sum(fixed > -0.001 & fixed <=-0.00001);
	catn(stats_P3 + "," + fixed_VstrDel_muts + "," + fixed_strDel_muts + "," + fixed_modDel_muts + "," + fixed_wkDel_muts + "\n");
	if(p3.individuals.size() < 2){
			sim.simulationFinished();
			cat("The population has gone extinct");
	}
}
EOM
