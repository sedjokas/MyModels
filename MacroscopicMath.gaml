/***
* Name: tb1
* Author: selain Kasereka
* Description: Simulate TB in a large scale, case of DRC
* Tags: Tag1, Tag2, TagN
*/
model tb1

global {
	int S_people <- 500 min: 500 max:1000000; //Number of sc=usceptible individuals
	int I_people <- 100 min:1 max: 541583; // infectious people (active TB)
	int Le_people <- 5; // Number of early latent TB people
	int Lf_people <- 5; // Number of late latent TB people
	int T_people <- 2; // Number of active TB people transferred people to another hospital 
	int K_people <- 100; // Nimber of active TB people who interrupt their treatment;
	int R1_people <- 0; // Number of healthy people after treatment process and tested;
	int R2_people <- 0; // Number of spontaneously healthy people and tested
	float Gamma<-10.0; // rate of recruitment in compartment S 
	float mu1 <-0.003; //rate of mortality not related to TB infection 
	float mu2 <-0.003; //rate of mortality related to TB infection  
	float beta <- 0.01; //rate of transferred people in a other hospital
	float gamma <- 0.85; //rate of recovered after treatment process
	float sigma<-0.02; //rate of spontaneously recovered
	float alpha <-0.4; //rate of contact
	float lambda <-0.008; // rate of transmission
	float p <-0.5; // proportion of people
	float v<-0.5; // rate of infected who interrupt their treatment 
	float q<-0.3; // rate of fast progression to active TB 
	float r<- 0.03; //rate of re-infection from R1 and R2 to Le
	float r1<- 0.02; //rate of re-infection from I to Le
	float r2<- 0.02; //rate of re-infection from T to Lf
	float r3<- 0.02; //rate of re-infection from K to Le
	float g1<- 0.85; //rate of recovered after treatment process from Le to R1
	float g2<- 0.01; // rate of spontaneously recovered from Le to R2
	float k1<-0.85; // rate of recovered after treatment process from Lf to R1 
	float k2<-0.01; //rate of spontaneously recovered from Lf to R2
	float h<-0.5; //rate of progression of TB infection, from Le to Lf
	float w<-0.2; // rate of slow progression of TB infection, from Lf to I 
	
	
	//Total population
	int N<- S_people + I_people + Le_people + Lf_people + T_people + K_people + R1_people + R2_people;
	geometry shape <- square(50);
	init {
		//Create node agent for ODE System
		create math_TB number: 1 {
			S <- float(S_people);
			I <- float(I_people);
			Le <- float(Le_people);
			Lf <- float(Lf_people);
			T <- float(T_people);
			K <- float(K_people);
			R1<- float(R1_people);
			R2 <- float(R2_people);
		}
	}

}


//Species node agent that will represent the ODE system
species math_TB {
	float t;    
	float S <- N-I-Le-Lf-T-K; 
	float I;
	float Le;
	float Lf;
	float T;
	float K ;
	float R1 ;
	float R2 ;

	equation SIR { 
		    diff(S,t) = (Gamma - alpha*lambda*p *S * I / N - alpha*lambda*(1-p)*S*I/N - mu1*S);
		    diff(Le,t) = (alpha*lambda*p *S * I/N + alpha*lambda*r*I*(R1 + R2) + r1*I + r2*T + r3*K - (mu1 + h + q+ g1+ g2)*Le);
			diff(Lf,t) = (h*Le - (mu1 + w + k1 + k2)*Lf);
			diff(I,t)= (w*Lf + q*Le - (r1 + gamma + beta + sigma + v + mu1 + mu2)*I + alpha*lambda*R1*I + alpha*lambda*R2*I + alpha*lambda*(1-p)*S*I/N);
			diff(R1,t)= (g1*Le + k1*Lf + gamma*I - alpha*lambda*r*R1*I - mu1*R1 -alpha*lambda*R1*I);
			diff(R2,t)= (sigma*I + k2*Lf + g2*Le - alpha*lambda*R2*I/N - alpha*lambda*r*R2*I - mu1*R2);
			diff(T,t)= (beta*I - (mu1 + mu2 + r2)*T);
			diff(K,t)= (v*I - (mu1 + mu2 + r3)*K);
			}
				
    	reflex solving {solve SIR method: rk4 step: 0.7;}
}

experiment Simulation_Math_TB type: gui {
	//population
	parameter 'Number of Susceptible' type: int var: S_people category: "Initial population";
	parameter 'Number of Active TB' type: int var: I_people category: "Initial population";
	parameter 'Number early latent TB' type: int var: Le_people category: "Initial population";
	parameter 'Number late latent TB' type: int var: Lf_people category: "Initial population";
	parameter 'Number of transfered people' type: int var: T_people category: "Initial population";
	parameter 'Number of TB people who interrupt treatment' type: int var: K_people category: "Initial population";
	parameter 'Number of Recovered people after treatment' type: int var: R1_people category: "Initial population";
	parameter 'Number of spontaneously Recovered people' type: int var: R2_people category: "Initial population";
	// my parameters
	parameter 'Gamma/recrut' type: float var: Gamma <- 0.1 category: "Parameters";
	parameter 'Mortality linked to TB' type: float var: mu1 <- 0.01 category: "Parameters";
	/*parameter 'Mortality not linked to TB' type: float var: mu2 <- 0.01 category: "Parameters";
	parameter 'Beta (S->I)' type: float var: beta <- 0.1 category: "Parameters";
	parameter 'Gamma (I->R)' type: float var: gamma <- 0.01 category: "Parameters";

	float Gamma<-0.0; // rate of recruitment in compartment S 
	float mu1 <-0.0; //rate of mortality not related to TB infection 
	float mu2 <-0.0; //rate of mortality related to TB infection  
	float beta <- 0.1; //rate of transferred people in a other hospital
	float gamma <- 0.0; //rate of recovered after treatment process
	float sigma<-0.0; //rate of spontaneously recovered
	float alpha <-0.0; //rate of contact
	float lambda <-0.0; // rate of transmission
	float p <-0.0; // proportion of people
	float v <-0.0; // rate of infected who interrupt their treatment 
	float q <-0.0; // rate of fast progression to active TB 
	float r <- 0.0; //rate of re-infection from R1 and R2 to Le
	float r1 <- 0.0; //rate of re-infection from I to Le
	float r2 <- 0.0; //rate of re-infection from T to Lf
	float r3 <- 0.0; //rate of re-infection from K to Le
	float g1 <- 0.0; //rate of recovered after treatment process from Le to R1
	float g2 <- 0.0; // rate of spontaneously recovered from Le to R2
	float k1 <-0.0; // rate of recovered after treatment process from Lf to R1 
	float k2 <-0.0; //rate of spontaneously recovered from Lf to R2
	float h <-0.0; //rate of progression of TB infection, from Le to Lf
	float w<-0.0; //  rate of slow progression of TB infection, from Lf to I */
	
	output {
		
		display TB1 { 
			chart "TB - STAT" type: series background: #white {
				data 'S' value: first(math_TB).S marker: true color: #green;
				data 'I' value: first(math_TB).I marker: true  color: #red;
				data 'Le' value: first(math_TB).Le marker: true color: #orange;
				data 'Lf' value: first(math_TB).Lf marker: true color: #yellow;
				data 'T' value: first(math_TB).T marker: true color: #gray;
				data 'K' value: first(math_TB).K marker: true color: #magenta;
				data 'R1' value: first(math_TB).R1 marker: true color: #blue;
				data 'R2' value: first(math_TB).R2 marker: true color:rgb(#77B5FE);
			}
		}
		
		display "datalist_pie_chart" type: java2D
		{
			chart "datalist_pie_chart" type: pie style: exploded
			{
				datalist legend: ["S", "I", "Le", "Lf", "T", "K", "R1", "R2"] value: [[first(math_TB).S], [first(math_TB).I],[first(math_TB).Le],[first(math_TB).Lf],[first(math_TB).T],
					[first(math_TB).K],[first(math_TB).R1],[first(math_TB).R2]
				] 
				color: [# green, # red, # orange, # yellow, # gray, # magenta, # blue, rgb(#77B5FE)];
			}

		}
				
		display TB2 { 
			chart "TB - STATCam" type: pie style: "3d"{
				data 'S' value: first(math_TB).S color: #green;
				data 'I' value: first(math_TB).I color: #red;
				data 'Le' value: first(math_TB).Le color: #orange;
				data 'Lf' value: first(math_TB).Lf color: #yellow;
				data 'T' value: first(math_TB).T color: #gray;
				data 'K' value: first(math_TB).K color: #magenta;
				data 'R1' value: first(math_TB).R1 color: #blue;
				data 'R2' value: first(math_TB).R2 color:rgb(#77B5FE);
			}
		}
	}

}


