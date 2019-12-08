/*****
 *  Modele-Hybride TB
 * 	Cadre: PhD
 *  Author: Selain Kasereka
 *  E-mail: selain.kasereka@unikin.ac.cd
 *  Supervision: Emile-Franc Doungmo Goufo and Ho Tuong Vinh
 *  Description: A hybrid meta-populationnal model that couple ABM and EBM, the dynamics of the disease is managed by 
    the differential equations solved by the rk4 method in each city. The bihavior of all andividuals is 
    considered separetly and managed by ABM. 
 *  -----------------------------------------------------------
 *  N: Total population 
 *  iInit: initial infected individual
 *  stepI: discrétisation step for the resolution of the ODE
 *****/
 
model HybridTB

global {
	/** Insert the global definitions, variables and actions here */
    
//My parameter
	float Gamma ; // Rate of recruitment in compartment S 
	float mu1 ; //rate of mortality not related to TB infection 
	float mu2 ; //rate of mortality related to TB infection  
	float beta ; //rate of transferred people in a other hospital
	float gamma ; //rate of recovered after treatment process
	float sigma ; //rate of spontaneously recovered
	float alpha ; //rate of contact
	float lambda ; // rate of transmission
	float p ; // proportion of people
	float v ; // rate of infected who interrupt their treatment 
	float q ; // rate of progression to active TB (Le ->I)
	float r ; //rate of re-infection from R1 and R2 to Le
	float r1 ; //rate of re-infection from I to Le
	float r2 ; //rate of re-infection from T to Lf
	float r3 ; //rate of re-infection from K to Le
	float g1 ; //rate of recovered after treatment process from Le to R1
	float g2 ; // rate of spontaneously recovered from Le to R2
	float k1 ; // rate of recovered after treatment process from Lf to R1 
	float k2 ; //rate of spontaneously recovered from Lf to R2
	float h ; //rate of progression of TB infection, from Le to Lf
	float w ; // rate of slow progression of TB infection, from Lf to I 
	
//My population
	int people <- 10000; //min: 1 max: 100000; //Population
	float step <- 0.01;//min: 0.0001 max: 0.1; //Step
	int iInit <- 1;// min: 1 max: 1000; //Initial number of infectious people
	int NbVille <- 6;// min: 1 max: 50; // Number of city
	float rate_S <- 0.01;// min: 0.01 max: 100.0; //Rate of creation of Susceptible people
	float rate_I <- 0.00001;// min: 0.01 max: 100.0; //rate of creation of infectious people
	float rate_Le <- 0.00001;// min: 0.01 max: 100.0; // rate of creation of recovered people
	float rate_Lf <- 0.00001;// min: 0.01 max: 100.0; //Rate of creation of Susceptible people
	float rate_T <- 0.01;// min: 0.01 max: 100.0; //rate of creation of infectious people
	float rate_K <- 0.01;// min: 0.01 max: 100.0; // rate of creation of recovered people
	float rate_R1 <- 0.01;// min: 0.01 max: 100.0; //Rate of creation of Susceptible people
	float rate_R2 <- 0.01;// min: 0.01 max: 100.0; //rate of creation of infectious people
	int idVille <- 0;
	int idville_act <- 0 ;//min:1 max:100;
	float max_I <- float(iInit);
	int k <- 0;
	int temps_Imax <- cycle;
	
	geometry shape <- square(50);
	
	init{
				
		let i<-0;
		let THETA<-360/NbVille;	
		let width<-50;
		let height<-50;
		
		create ville number: NbVille {
			compte_ville <- idVille;
			idVille <- idVille + 1;
			
			 if (idVille =1)
			 {
			 	I <- iInit;
			 	S <- people / NbVille - iInit;
			 }
//write "ville numero " + idVille + "  il ya infected " + I;
			
			let r1<-(width - size_vil*2)/2;
					
				location <-point([width/2 + r1 * cos(THETA * i), height/2 + r1 * sin(THETA *i )]);
				i<-i+1;	
						
			}
		
		create hub number:1{
			location<- point([width/2,height/2]);
		}
		
		create statistique number:1{
			
		}
	save "cycle,Infected" to:"fileoutM2.txt";

}

// Ecrire le resultat de la simulation dans le fichier "fileout.txt" (voir repertoire models)
 
 /*reflex stop when: (sum(ville collect each.I) < 1) {
		write " I max "+ max_I;
		write "Nombre de cycle a I max "+ temps_Imax;
		write "Durée de l'épidémie " + cycle + " cycles";
		write "------------------------------------------ ";
		
		do halt;
		
	} */

reflex stop_simulation when: (cycle = 700) {
					do pause ;
			} 
		
reflex save_data when: every(1#cycle){			 	
		
		//save the following text into the given text file. Note that each time the save statement is used, a new line is added at the end of the file.
		save [cycle,sum(ville collect each.S) ,sum(ville collect each.I) ,sum(ville collect each.Le) ,sum(ville collect each.Lf) ,sum(ville collect each.T) ,sum(ville collect each.K) ,sum(ville collect each.R1) ,sum(ville collect each.R2) ] to: "../results/dataAllcity.csv" type:"csv" rewrite: false ;
		save [cycle,first((ville as list) where (each.compte_ville=0)).S,first((ville as list) where (each.compte_ville=0)).I,first((ville as list) where (each.compte_ville=0)).Le,first((ville as list) where (each.compte_ville=0)).Lf,first((ville as list) where (each.compte_ville=0)).T,first((ville as list) where (each.compte_ville=0)).K,first((ville as list) where (each.compte_ville=0)).R1,first((ville as list) where (each.compte_ville=0)).R2] to: "../results/dataCity1.csv" type:"csv" rewrite: false;
		save [cycle,first((ville as list) where (each.compte_ville=1)).S,first((ville as list) where (each.compte_ville=1)).I,first((ville as list) where (each.compte_ville=1)).Le,first((ville as list) where (each.compte_ville=1)).Lf,first((ville as list) where (each.compte_ville=1)).T,first((ville as list) where (each.compte_ville=1)).K,first((ville as list) where (each.compte_ville=1)).R1,first((ville as list) where (each.compte_ville=1)).R2] to: "../results/dataCity2.csv" type:"csv" rewrite: false;
		save [cycle,first((ville as list) where (each.compte_ville=2)).S,first((ville as list) where (each.compte_ville=2)).I,first((ville as list) where (each.compte_ville=2)).Le,first((ville as list) where (each.compte_ville=2)).Lf,first((ville as list) where (each.compte_ville=2)).T,first((ville as list) where (each.compte_ville=2)).K,first((ville as list) where (each.compte_ville=2)).R1,first((ville as list) where (each.compte_ville=2)).R2] to: "../results/dataCity3.csv" type:"csv" rewrite: false;
		save [cycle,first((ville as list) where (each.compte_ville=3)).S,first((ville as list) where (each.compte_ville=3)).I,first((ville as list) where (each.compte_ville=3)).Le,first((ville as list) where (each.compte_ville=3)).Lf,first((ville as list) where (each.compte_ville=3)).T,first((ville as list) where (each.compte_ville=3)).K,first((ville as list) where (each.compte_ville=3)).R1,first((ville as list) where (each.compte_ville=3)).R2] to: "../results/dataCity4.csv" type:"csv" rewrite: false;
		save [cycle,first((ville as list) where (each.compte_ville=4)).S,first((ville as list) where (each.compte_ville=4)).I,first((ville as list) where (each.compte_ville=4)).Le,first((ville as list) where (each.compte_ville=4)).Lf,first((ville as list) where (each.compte_ville=4)).T,first((ville as list) where (each.compte_ville=4)).K,first((ville as list) where (each.compte_ville=4)).R1,first((ville as list) where (each.compte_ville=4)).R2] to: "../results/dataCity5.csv" type:"csv" rewrite: false;
		save [cycle,first((ville as list) where (each.compte_ville=5)).S,first((ville as list) where (each.compte_ville=5)).I,first((ville as list) where (each.compte_ville=5)).Le,first((ville as list) where (each.compte_ville=5)).Lf,first((ville as list) where (each.compte_ville=5)).T,first((ville as list) where (each.compte_ville=5)).K,first((ville as list) where (each.compte_ville=5)).R1,first((ville as list) where (each.compte_ville=5)).R2] to: "../results/dataCity6.csv" type:"csv" rewrite: false;	    	     
	} 
}

species ville {
		
	rgb color_vil <-rgb('yellow'); 
	int size_vil <-6;
	int N <- round(people/NbVille);
	int Ninit<-N;
	int NSortis<-0;
	int compte_ville;

  	int dS<-0;
  
  	float t;    
	float S <- N-I-Le-Lf-T-K; 
	float I;
	float Le;
	float Lf;
	float T;
	float K ;
	float R1 ;
	float R2 ;
			
/* 	equation SLeLfITKR1R2 { 
		    diff(S,t) = (Gamma - alpha*lambda*p *S * I/N - alpha*lambda*(1-p)*S*I/N - mu1*S);
		    diff(Le,t) = (alpha*lambda*p *S * I/N + alpha*lambda*r*I*(R1 + R2)/N + r1*I + r2*T + r3*K - (mu1 + h + q+ g1+ g2)*Le);
			diff(Lf,t) = (h*Le - (mu1 + w + k1 + k2)*Lf);
			diff(I,t) = (w*Lf + q*Le - (r1 + gamma + beta + sigma + v + mu1 + mu2)*I + alpha*lambda*R1*I/N + alpha*lambda*R2*I/N + alpha*lambda*(1-p)*S*I/N);
			diff(R1,t) = (g1*Le + k1*Lf + gamma*I - alpha*lambda*r*R1*I/N - mu1*R1 -alpha*lambda*R1*I/N);
			diff(R2,t) = (sigma*I + k2*Lf + g2*Le - alpha*lambda*R2*I/N - alpha*lambda*r*R2*I/N - mu1*R2);
			diff(T,t) = (beta*I - (mu1 + mu2 + r2)*T);
			diff(K,t) = (v*I - (mu1 + mu2 + r3)*K);
			}			
	 */	 
	 			
		equation SLeLfITKR1R2 { 
		    diff(S,t) = (Gamma*N - alpha*lambda*p *S * I - alpha*lambda*(1-p)*S*I - mu1*S);
		    diff(Le,t) = (alpha*lambda*p *S * I + alpha*lambda*r*I*(R1 + R2) + r1*I + r2*T + r3*K - (mu1 + h + q+ g1+ g2)*Le);
			diff(Lf,t) = (h*Le - (mu1 + w + k1 + k2)*Lf);
			diff(I,t) = (w*Lf + q*Le - (r1 + gamma + beta + sigma + v + mu1 + mu2)*I + alpha*lambda*R1*I + alpha*lambda*R2*I + alpha*lambda*(1-p)*S*I);
			diff(R1,t) = (g1*Le + k1*Lf + gamma*I - alpha*lambda*r*R1*I - mu1*R1 -alpha*lambda*R1*I);
			diff(R2,t) = (sigma*I + k2*Lf + g2*Le - alpha*lambda*R2*I - alpha*lambda*r*R2*I - mu1*R2);
			diff(T,t) = (beta*I - (mu1 + mu2 + r2)*T);
			diff(K,t) = (v*I - (mu1 + mu2 + r3)*K);
			}
	
			
	    reflex solving  {
    					solve SLeLfITKR1R2 method:rk4 step: step;
    					}	
			
	aspect voir_ville {
		
		draw circle(size_vil) color: color_vil;
			color_vil <- rgb([I,1.0,1.0]);			
			
			}
			
 int sain<-0;
 int infected <-0;
 int immune<-0;

 int cpt_creer_homme<-1;
 int cpt_voyage_homme<-50;

 int faire_action<-1; // 1 pour oui et 0 pour non
 int dcmpte<-cpt_creer_homme;
	
 reflex cpteur {

  dcmpte<-dcmpte -1;
 	if(dcmpte=0 and faire_action=1){
 		do creer_homme_in_ville;
 		faire_action<-0;
 		dcmpte <-cpt_voyage_homme;
 		
 	}
 	else if(dcmpte=0 and faire_action=0)
 	{
 		dcmpte<-cpt_creer_homme;
 		faire_action<-1;
 	}
 }
  
 
action creer_homme_in_ville {
	if (NSortis<=Ninit) {
 			if compte_ville=1{

 	}
 
 			
//Create infected people and select the reate of travel people based on theier status (S)
 			
 			sain<-floor(S*rate_S);				
 			if (Ninit-NSortis<sain){
 				sain<-Ninit-NSortis;
 			}
 			if (Ninit-NSortis>=sain){
 					create homme number:sain  {
			     	 	location <- myself.location;
			     	 	villeDepart<-myself;
			       	   	is_susceptible <- true;
			        	is_infected <-  false;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <-  false; 
			            is_immune2 <-  false; 
			            is_lost <-  false;
			            is_trans <-  false;
			            color <-  rgb('green');
			            depart<-true;
			           
 						}					
 					dS<-dS+sain;
 					NSortis<-NSortis+floor(float(sain));
			}
			
//Create infected people and select the reate of travel people based on theier status (I) 
 
 			infected<-floor(I*rate_I);
 			// write "inf3    " + I
 			if (Ninit-NSortis<infected){infected<-Ninit-NSortis; }
 				if (Ninit-NSortis>=infected){
 					create homme number: infected {
			     	 	set location <- myself.location;
			     	 	villeDepart<-myself;
			       	   	is_susceptible <- false;
			        	is_infected <-  true;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <-  false; 
			            is_immune2 <-  false; 
			            is_lost <-  false;
			            is_trans <-  false; 
			            color <-  rgb('red');
			            depart<-true;
			            
            		 }
            		 //set dS<-infected;
            		    NSortis<-NSortis+floor(float(infected));
				}

 //Create recovered after treatment people and select the reate of travel people based on theier status (R1)
  
			let immune1<-floor(R1*rate_R1);
			//write "immu3    " + R;
			if (Ninit-NSortis<immune1){immune1<-Ninit-NSortis; }
				if (Ninit-NSortis>=immune1){
 					create homme number: immune1 {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <-  false;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <-  true; 
			            is_immune2 <-  false; 
			            is_lost <-  false;
			            is_trans <-  false; 
			            color <- #blue;
			            depart<-true;
			            
	 					}
 				//set dR<-immune;
 				NSortis<-NSortis+floor(float(immune1));
				}
			// write "les retires qui partent  " + floor(float(immune));

 //Create recovered spontaneously people and select the reate of travel people based on rate_R2 
 			
			let immune2<-floor(R2*rate_R2);
			//write "immu3    " + R;
			if (Ninit-NSortis<immune2){
				immune2<-Ninit-NSortis;
			}
				if (Ninit-NSortis>=immune2){
 					create homme number: immune2 {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <- false;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <- false; 
			            is_immune2 <- true; 
			            is_lost <- false;
			            is_trans <- false; 
			            color <-  rgb(#77B5FE);
			            depart<-true;
			            
	 					}
 				
 				NSortis<-NSortis+floor(float(immune2));
				}
			// write "les retires qui partent  " + floor(float(immune2));
			
 //Create Latent Late people and select the reate of travel people based on rate_LK  	
 	
 		let transfert<-floor(T*rate_T);
			
			if (Ninit-NSortis<transfert){transfert<-Ninit-NSortis; }
				if (Ninit-NSortis>=transfert){
 					create homme number: transfert {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <- false;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <- false; 
			            is_immune2 <- false; 
			            is_lost <- false;
			            is_trans <- true; 
			            color <- #gray;
			            depart<-true;
			            
	 					}
 				
 				NSortis<-NSortis+floor(float(transfert));
				}

			
 //Create Latent Late people and select the reate of travel people based on their status (K) 
 			
 		let lost<-floor(K*rate_K);
			
			if (Ninit-NSortis<lost){set lost<-Ninit-NSortis; }
				if (Ninit-NSortis>=lost){
 					create homme number: lost {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <- false;
			            is_Le <-  false; 
			            is_Lf <- false;
			            is_immune1 <- false; 
			            is_immune2 <- false; 
			            is_lost <- true;
			            is_trans <- false; 
			            color <- #magenta;
			            depart<-true;
			            
	 					}
 				
 				NSortis<-NSortis+floor(float(lost));
				}
				
 //Create Latent Late people and select the reate of travel people based on their status (Le)  	
 	let LatentEarly<-floor(Le*rate_Le);
			
			if (Ninit-NSortis<LatentEarly){LatentEarly<-Ninit-NSortis; }
				if (Ninit-NSortis>=LatentEarly){
 					create homme number: LatentEarly {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <- false;
			            is_Le <-  true; 
			            is_Lf <- false;
			            is_immune1 <- false; 
			            is_immune2 <- false; 
			            is_lost <- false;
			            is_trans <- false; 
			            color <- #orange;
			            depart<-true;
			            
	 					}
 				
 				NSortis<-NSortis+floor(float(LatentEarly));
				}
				
 //Create Latent Late people and select the reate of travel people based on their status (Lf) 
 	
 	let LatentLate<-floor(Lf*rate_Lf);
			
			if (Ninit-NSortis<LatentLate)
			{
				LatentLate<-Ninit-NSortis;
			}
				if (Ninit-NSortis>=LatentLate){
 					create homme number: LatentLate {
			     	 	location<-myself.location;
			     	 	villeDepart<-myself;
			     	  	is_susceptible <- false;
			        	is_infected <- false;
			            is_Le <- false; 
			            is_Lf <- true;
			            is_immune1 <- false; 
			            is_immune2 <- false; 
			            is_lost <- false;
			            is_trans <- false; 
			            color <- #yellow;
			            depart<-true;
			            
	 					}
 				
 				NSortis<-NSortis+floor(float(LatentLate));
				}
 	
 	
 		} 
	}      


}

species hub{
	rgb color_hub<-#black;
	int size_hub<-1;
	
	aspect voir_hub {
			draw circle(size_hub ) color: color_hub;
			
			}
	
}
		
species homme skills:[moving] {
		
	    float size_hom <- 0.3; 
		bool is_susceptible <- true;
		bool is_infected <-  false;
		bool is_Le <-  false; 
		bool is_Lf <- false;
		bool is_immune1 <-  false; 
		bool is_immune2 <-  false; 
		bool is_lost <-  false;
		bool is_trans <-  false;
        ville villeDepart;
        rgb color;
        ville villeDest<-one_of(list(ville)-self);
        point dest <-point(villeDest);
        point passage<-point(one_of(hub));
        point objectif<-dest;
	    int sain <-0;
	    //int infected <-0;
	    //int immune <-0;
	    bool depart<-false;
	    int decompter<-0;
	    
		aspect voir_hom {
			draw circle(size_hom) color: color; 
			}	

	reflex voyager {	
		
		if(depart=true){
			if(decompter=0){
					decompter<-1;
					if (is_susceptible){
						villeDepart.S<-villeDepart.S-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dS<-villeDepart.dS+1;
					}
					else if(is_infected){
						villeDepart.I<-villeDepart.I-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else if(is_Le){
						villeDepart.Le <-villeDepart.Le-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else if(is_Lf){
						villeDepart.Lf <-villeDepart.Lf-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else if(is_lost){
						villeDepart.K <-villeDepart.K-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else if(is_trans){
						villeDepart.T <-villeDepart.T-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else if(is_immune1){
						villeDepart.R1 <-villeDepart.R1-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dI<-villeDepart.dI+1;
					}
					else{
						 villeDepart.R2<-villeDepart.R2-1;
						villeDepart.N<-villeDepart.N-1;
						//set villeDepart.dR<-villeDepart.dR+1;
					}
			}
			if (location  != passage and passage!=nil)
			{
				do goto target: passage speed:50;	
			}
			else if(location  = passage and passage!=nil){
				set passage<-nil;
				do goto target: objectif speed:50;
			}
			else if(location  != objectif and passage=nil and objectif!=nil){
				do goto target: objectif speed:50;
			}
			else if(location = objectif and passage=nil and objectif!=nil){
				if (is_susceptible){
					villeDest.S<-villeDest.S+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				else if(is_infected){
					villeDest.I<-villeDest.I+1;
					villeDest.N<-villeDest.N+1;
					//set villeDest.aI<-villeDest.aI-1;
				}
				if (is_Le){
					villeDest.Le<-villeDest.Le+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				if (is_Lf){
					villeDest.Lf<-villeDest.Lf+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				if (is_lost){
					villeDest.T<-villeDest.T+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				if (is_trans){
					villeDest.T<-villeDest.T+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				if (is_immune1){
					villeDest.R1<-villeDest.R1+1;
					villeDest.N<-villeDest.N+1;	
					//set villeDest.aS<-villeDest.aS-1;				
				}
				if (is_immune2){
					villeDest.R2<-villeDest.R2+1;
					villeDest.N<-villeDest.N+1;
					//set villeDest.aR<-villeDest.aR-1;
				}
				
				set objectif<-nil;
				//do die;
			}
	
		}	
	}		

}

species statistique {
	float som_S<- float(people - iInit);
	float som_I<-float(iInit);
	float som_Le<-0.0;
	float som_Lf<-0.0;
	float som_T<-0.0;
	float som_K<-0.0;
	float som_R1<-0.0;
	float som_R2<-0.0;
	float som_N <-som_S + som_I + som_Le + som_Lf + som_R1 + som_R2 + som_K + som_T;

reflex ecrire_file{
 	
 	//save ""+cycle+","+ som_I  type:"csv" to:"fileoutM2.txt" rewrite:false;
 
 }
	reflex recalcul_SI {
		som_S <- sum(ville collect each.S) + length(homme where each.is_susceptible); 
		som_I <- sum(ville collect each.I) + length(homme where each.is_infected);
		som_Le <- sum(ville collect each.Le) + length(homme where each.is_Le);
		som_Lf <- sum(ville collect each.Lf) + length(homme where each.is_Lf); 
		som_T <- sum(ville collect each.T) + length(homme where each.is_trans);
		som_K <- sum(ville collect each.K) + length(homme where each.is_lost);
		som_R1 <- sum(ville collect each.R1) + length(homme where each.is_immune1); 
		som_R2 <- sum(ville collect each.R2) + length(homme where each.is_immune2);
		som_N<- som_S + som_I + som_Le + som_Lf + som_R1 + som_R2 + som_K + som_T; 
		
		//write "toto " + length(homme where each.is_susceptible);
		
		if (som_I>max_I){
			set max_I<-som_I;
			set temps_Imax<-cycle;
			
			}
		}
	}


experiment HybridTB type: gui {

	// my population
	parameter 'All people' type: int var: people category: 'Equation SIR'; 
	parameter 'Nomber of city' type: int var: NbVille category: 'Ville';
	parameter 'City ID' type: int var: idville_act category: 'Ville'; 
	
	// my parameters	   
	parameter 'Recrutment: Gamma' type: float var: Gamma <- 0.0300 category: "Parameters"; //Recrutment
	parameter 'Natural Mortality: mu1' type: float var: mu1 <- 0.0222 category: "Parameters";//rate of mortality related to TB infection  
	parameter 'Mortality linked to TB: mu2' type: float var: mu2 <- 0.040  category: "Parameters"; //rate of mortality not related to TB infection 0.0003
	parameter 'Rate of transfer: beta (I->T)' type: float var:beta <- 0.010 category: "Parameters"; //rate of transferred people in a other hospital
	parameter 'Rate of recovered after treatment: gamma (I->R1)' type: float var:gamma <- 0.84 category: "Parameters"; //rate of recovered after treatment process
	parameter 'Rate of spontaneously recovered: sigma (I->R2)' type: float var: sigma<-0.25 category: "Parameters"; //rate of spontaneously recovered
	parameter 'Rate of contact' type: float var:alpha <-0.10 category: "Parameters"; //rate of contact
	parameter 'Rate of transmission: lambda' type: float var: lambda <- 0.10 category: "Parameters"; // rate of transmission
	parameter 'Proportion of people: p (S-> Le & I)' type: float var: p <- 0.95  category: "Parameters"; // proportion of people
	parameter 'Rate of progression to ATB: q (Le->I)' type: float var: q <-0.129  category: "Parameters"; // Rate of progression to ATB: q (Le->I)
	parameter 'Rate of progression to ATB: w (Lf->I)' type: float var: w <-0.075  category: "Parameters"; // Rate of progression to ATB: w (Lf->I)
	parameter 'Rate of progression to Latent late: h (Le->Lf)' type: float var:h <- 0.821  category: "Parameters"; //rate of re-infection from R1 and R2 to Le
	parameter 'Rate of re-infection: r1 (I->Le)' type: float var: r1 <- 0.63  category: "Parameters"; //rate of re-infection from I to Le
	parameter 'Rate of re-infection: r2 (T->Le)' type: float var: r2 <- 0.63  category: "Parameters"; //rate of re-infection from T to Le
	parameter 'Rate of re-infection: r3 (K->Le)' type: float var: r3 <- 0.63  category: "Parameters"; //rate of re-infection from K to Le
	parameter 'Recovered after treatment: g1 (Le->R1)' type: float var: g1 <- 0.84  category: "Parameters"; //rate of recovered after treatment process from Le to R1
	parameter 'Recovered spontaneously: g2 (Le->R2)' type: float var: g2 <- 0.25  category: "Parameters"; // rate of spontaneously recovered from Le to R2
	parameter 'Recovered after treatment: k1 (Lf->R1)' type: float var: k1 <- 0.84 category: "Parameters"; //rate of recovered after treatment process from Lf to R1
	parameter 'Recovered spontaneously: k2(Lf->R2)' type: float var: k2 <- 0.25  category: "Parameters"; // rate of spontaneously recovered from Lf to R2
	parameter 'Rate of reactivation: r (R1 & R2 -> Le)' type: float var: r <- 0.030  category: "Parameters"; // Rate of reactivation: r (R1 & R2 -> Le)
	parameter'Rate of treatment interuption: v (I->K)' type: float var: v <- 0.03  category: "Parameters"; //Rate of treatment interuption v (I->K)

 // Parametre a manupuler par les utilisateurs

	parameter 'Initial Infected people' var: iInit category: 'Equation SIR'; 
	parameter 'Step' var: step category: 'Equation SIR';
			 
	output { 

// Graphique pour toutes les villes a la fois, donc toute la population
				
      display voir refresh_every: 1 {
     // 	 grid model_grid lines: #black;
      	 species homme aspect:voir_hom;
      	 species ville aspect:voir_ville;
      	 species hub aspect:voir_hub;
		 
		 		
		}

		display All_City_Serie refresh_every: 1 {
			chart "" type: series background: rgb('white') {
				data 'S' value: sum(ville collect each.S) color: #green ;				
				data 'I' value: sum(ville collect each.I) color: #red;
				data 'Le' value: sum(ville collect each.Le) color: #orange;
				data 'Lf' value: sum(ville collect each.Lf) color: #yellow ;				
				data 'T' value: sum(ville collect each.T) color: #gray ;
				data 'K' value: sum(ville collect each.K) color: #magenta ;
				data 'R1' value: sum(ville collect each.R1) color: #blue ;				
				data 'R2' value: sum(ville collect each.R2) color: rgb(#77B5FE) ;
			}
		}
		
		display All_City_Camembert refresh_every: 1 {
			chart "CITY 1" type: pie background: rgb('white') {
				data 'S' value: sum(ville collect each.S) color: #green ;				
				data 'I' value: sum(ville collect each.I) color: #red;
				data 'Le' value: sum(ville collect each.Le) color: #orange;
				data 'Lf' value: sum(ville collect each.Lf) color: #yellow ;				
				data 'T' value: sum(ville collect each.T) color: #gray ;
				data 'K' value: sum(ville collect each.K) color: #magenta ;
				data 'R1' value: sum(ville collect each.R1) color: #blue ;				
				data 'R2' value: sum(ville collect each.R2) color: rgb(#77B5FE) ;
			}
		}
		
	/*display All_city_camambert refresh_every: 1 {
			chart "SIR_MATHS_MODEL2" type: pie {
				data 'S' value: first(statistique).som_S color: #green ;					
				data 'I' value: first(statistique).som_I color: #red;
				data 'Le' value: first(statistique).som_Le  color: #orange ;
				data 'Lf' value: first(statistique).som_Lf color: #yellow ;					
				data 'T' value: first(statistique).som_T color: #gray ;
				data 'K' value: first(statistique).som_K color: #magenta ;					
				data 'R1' value: first(statistique).som_R1 color: #blue;
				data 'R2' value: first(statistique).som_R2  color: rgb(#77B5FE)  ;
			}
		}*/
			display City_1 refresh_every: 1 {
			chart "CITY 1" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=0)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=0)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=0)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=0)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=0)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=0)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=0)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=0)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
			display City_2 refresh_every: 1 {
			chart "SIR_MATHS_" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=1)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=1)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=1)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=1)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=1)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=1)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=1)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=1)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
			
	display City_3 refresh_every: 1 {
			chart "SIR_MATHS_" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=2)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=2)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=2)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=2)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=2)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=2)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=2)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=2)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
		
			display City_4 refresh_every: 1 {
			chart "SIR_MATHS_" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=3)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=3)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=3)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=3)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=3)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=3)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=3)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=3)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
		display City_5 refresh_every: 1 {
			chart "SIR_MATHS_" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=4)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=4)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=4)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=4)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=4)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=4)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=4)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=4)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
	
	display City_6 refresh_every: 1 {
			chart "SIR_MATHS_" type: series background: rgb('white')  {
				data 'S' value: first((ville as list) where (each.compte_ville=5)).S color: #green ;				
				data 'I' value: first((ville as list) where (each.compte_ville=5)).I color:  #red;
				data 'Le' value: first((ville as list) where (each.compte_ville=5)).Le color: #orange ;
				data 'Lf' value:first((ville as list) where (each.compte_ville=5)).Lf color:  #yellow  ; 
				data 'T' value: first((ville as list) where (each.compte_ville=5)).T color: #gray ; 
				data 'K' value: first((ville as list) where (each.compte_ville=5)).K color:  #magenta  ;
				data 'R1' value: first((ville as list) where (each.compte_ville=5)).R1 color: #blue ;
				data 'R2' value: first((ville as list) where (each.compte_ville=5)).R2 color: rgb(#77B5FE)  ;
				
			}
		}
	}
 
 }

