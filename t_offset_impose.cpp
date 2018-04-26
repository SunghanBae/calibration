int t_offset_impose(TString fin){

		Long64_t t,pret,ext,extpret,delt,delext,evt,extTimestamp[10000],extTprev,old_ext,extT0;
		Long64_t t_offset1,t_offset2,t_offset,extT0prev,jump,ext_keep,jumpleap;

		time_t tm;
		int time0,multinum;
		Short_t rangeType,infoCode;
		Int_t adcData,feeNo,chNo,quality,dssdNo,stripNo;
		double offset,gain,adcData_off,thresh;
		bool endflag,addflag,nextflag,rewind_flag,keepflag,leapflag;
		int i,j,asicNo,asicNo0,inext,i0,dssd0;

		TFile* f = new TFile(fin);
		TTree* tr = (TTree*)f->Get("tree");
		tr->SetBranchAddress("evt",&evt);
		tr->SetBranchAddress("timestamp",&t);	
		tr->SetBranchAddress("extTimestamp",&ext);	
		tr->SetBranchAddress("feeNo",&feeNo);	
		tr->SetBranchAddress("chNo",&chNo);	
		tr->SetBranchAddress("adcData",&adcData);
		tr->SetBranchAddress("dssdNo",&dssdNo);
		tr->SetBranchAddress("stripNo",&stripNo);
		tr->SetBranchAddress("rangeType",&rangeType);	
		tr->SetBranchAddress("infoCode",&infoCode);	

		TFile* cf= new TFile("$AIDADIR/calibration/aida_calibration_20180306_r59.root");
		TTree* ct = (TTree*)cf->Get("tree");
		ct->SetBranchAddress("gain",&gain);
		ct->SetBranchAddress("offset",&offset);
		ct->SetBranchAddress("quality",&quality);
		ct->SetBranchAddress("thresh",&thresh);


		int l = fin.Length();
		fin.Remove(l-5);
		TString fo= fin;
		//	fo+="_offset.root";							// no change in quality
		//	fo+="_offset_v2.root";						// change quality as 0 if adc <1500 for rangeType==0
		//	fo+="_offset_v3.root";						// change quality as 0 if adc >1500 for rangeType==0
//		fo+="_tsj_offset_v4.root";						// offset for noisy strips (especially, those in the edge.)
		//fo+="_tsj_offset_v5.root";						// offset for noisy strips and gain rough matching (especially, those in the edge.)
//		fo+="_tsj_offset_v6.root";						// offset and gain rough matching using R59
		fo+="_tsj_offset_v6_ASIC.root";						// offset and gain rough matching using R59, 200timestamp shift based on ASIC distinction

		TFile* nf = new TFile(fo,"RECREATE");
		TTree* newt = new TTree("tree","");
		newt->Branch("evt",&evt,"evt/L");
		newt->Branch("timestamp",&t,"timestamp/L");
		newt->Branch("extTimestamp",&ext,"extTimestamp/L");
		newt->Branch("old_extTimestamp",&old_ext,"old_extTimestamp/L");
		newt->Branch("feeNo",&feeNo,"feeNo/I");
		newt->Branch("chNo",&chNo,"chNo/I");
		newt->Branch("dssdNo",&dssdNo,"dssdNo/I");
		newt->Branch("stripNo",&stripNo,"stripNo/I");
		newt->Branch("adcData",&adcData,"adcData/I");
		newt->Branch("adcData_off",&adcData_off,"adcData_off/D");
		newt->Branch("rangeType",&rangeType,"rangeType/S");
		newt->Branch("infoCode",&infoCode,"infoCode/S");
		newt->Branch("quality",&quality,"quality/I");
		newt->Branch("gain",&gain,"gain/D");
		newt->Branch("offset",&offset,"offset/D");

		int n = tr->GetEntriesFast();

		std::vector<Long64_t> tmap,checker;
		std::vector<pair<Long64_t, Long64_t>> tmap2;
		
		bool dropflag = false;
		bool first_offset_flag = false;
				
		time0=time(&tm);	
		std::cout<<endl<<"AIDA extTimestamp - timestamp offset calculation"<<endl;
		
		tr->GetEntry(0);
		pret=ext;
		
		Long64_t temp_offset;
//		TTree* temptree=new TTree("temp","");
//		temptree->Branch("temp_offset",&temp_offset,"temp_offset/L");
		std::map<Long64_t, Long64_t> offset_map;
		std::map<Long64_t, Long64_t>::iterator ioffset_map;

		int offset_contain[100]={0,};
		int entry_contain[100]={0,};
		
		j=0;
		int starti;
		
		for(i=0;i<1000000;i++){

				tr->GetEntry(i);

				if(ext==-9999){
					continue;
				}else{
					temp_offset = ext-t;
					
					ioffset_map=offset_map.find(temp_offset);

					if(ioffset_map==offset_map.end()){
						std::cout<<"New offset found"<<endl;
						offset_map.insert(std::pair<Long64_t, Long64_t> (temp_offset,j));
						entry_contain[j]=i;
						offset_contain[j]+=1;
						j++;
					}else{
						offset_contain[ioffset_map->second]+=1;
						if(offset_contain[ioffset_map->second] > 800000) {
								t_offset1 = ioffset_map->first;
								starti	  = entry_contain[ioffset_map->second];
								std::cout<<"The timestamp offset obtained : "<<temp_offset<<" with "<<offset_contain[ioffset_map->second]<<" events out of "<<i<<". starti = "<<starti<<endl;
								break;
						}
					}
					
				}
				
				if(i%10000==0) {std::cout<<"\r"<<i<<"'s entry done. Timespent "<<time(&tm)-time0<<"s";std::cout.flush();}
		}
		
		std::cout<<endl<<"AIDA big timestamp jump correction"<<endl;

		tr->GetEntry(starti);
		extpret=ext;
		pret=t;
		keepflag=false;
		dropflag=false;
		leapflag=false;

		jumpleap=0;
		jump=0;


		for(i=starti;i<n;i++){

				if(i%100000==0) {std::cout<<"\r"<<i<<"'s entry done. Timespent "<<time(&tm)-time0<<"s";std::cout.flush();}

				tr->GetEntry(i);

				if(dropflag&&leapflag){ cout<<endl<<"Two flag fired at the same time! at "<<i<<endl; return -1;}
				
				if(!dropflag && t-pret > 1e8 ){							//timestamp leap occurs

						leapflag=true;
						jumpleap=(Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));
						ext= t + t_offset1 + jumpleap;

				}
				else if(dropflag && t-pret> 1e8){

						if((t-pret) > (jump+1e8)){
								
							jumpleap= jump + (Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));						//leap during the drop
							ext= t+ t_offset1 + jumpleap;
							leapflag=true;
							dropflag=false;
							
//							cout<<ext<<" ext and entry "<<i<<" with jumpleap "<<jumpleap<<endl; return -1;
										
						}else if((t-pret) < (jump-1e8) ){																	//drop partially recovered
							
							jump= jump + (Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));
							ext= t+ t_offset1 + jump;
							leapflag=false;
							dropflag=true;
							
						}else {																							//timestamp is recovered from drop	
							ext= t+ t_offset1;
							leapflag=false;
							dropflag=false;
							jump=0;
							jumpleap=0;																			
						}

				}else if(!leapflag && t-pret < -1e8){ 																	//timestamp drop occurs
				
						jump=(Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));
						dropflag=true;
						ext= t + t_offset1 + jump;
		
				}else if(leapflag && t-pret < -1e8){ 		

						if((t-pret) > (jumpleap +1e8)){																	//leap partially recovered only
						
							jumpleap= jumpleap + (Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));
							ext= t+ t_offset1 + jumpleap;
							leapflag=true;
							dropflag=false;

						}else if((t-pret) < (jumpleap - 1e8)){															//drop occurs right after the leap
							
							jump= jumpleap + (Long64_t)(ceil(1.0*(pret-t)/pow(2,28))*pow(2,28));
							ext= t+ t_offset1 + jump;
						//	cout<<ext<<" ext and entry "<<i<<" with jump "<<jump<<endl; return -1;
							leapflag=false;
							dropflag=true;

						}else{																							//timestamp recovered from leap
								
								dropflag=false;
								leapflag=false;
								ext= t + t_offset1;
								jump=0;
								jumpleap=0;

						}

				}else{
						if(leapflag){
								if(jumpleap >-1e8) {cout<<endl<<"The leap happenned. But no suitable jumpleap exist at i "<<i<<endl; return -1;}
								else{ext= t + t_offset1 + jumpleap;}
						}else if(dropflag){
								if(jump < 1e8) {cout<<endl<<"The drop happenned. But no suitable jump exist at i "<<i<<endl; return -1;}
								else{ext= t + t_offset1 + jump;}
						}else{
								ext= t + t_offset1;
						}

				}

//				if(ext<10){ cout<<i<<"entry has weird ext"<<endl; return -1;}

				tmap.push_back(ext);
				pret=t;
				checker.push_back(-1);
		}


		std::cout<<endl<<"ASIC multiplexer timestamp correction"<<endl;

		int iprev;
		int r_i;
		r_i=0;
		i=starti;
		j=0;
		rewind_flag=false;
		while(i<n){
				endflag=false;
				addflag=true;
				nextflag=false;
				multinum=0;

				tr->GetEntry(i);
				ext=tmap[i-starti];

				extT0=ext;

				asicNo0=(feeNo-1)*4 + (chNo/16);

				if(chNo<0){ 
						asicNo0=-1;
						tmap2.push_back(pair<Long64_t,Long64_t> (i,ext));
						i++; 
				}

				dssd0=dssdNo;

				while(!endflag && i<n && asicNo0>-1){

						if(addflag){
								tmap2.push_back(pair<Long64_t,Long64_t> (i,extT0));
								extTprev=ext;
								iprev=i;
								checker[i-starti]=1;
								addflag=false;
						}

						i++;

						if(checker[i-starti]!=1){

								tr->GetEntry(i);
								asicNo=(feeNo-1)*4+ (chNo/16);
								ext=tmap[i-starti];

								if(chNo<0){
										tmap2.push_back(pair<Long64_t,Long64_t>(i,extT0));
										checker[i]=1;
										continue;
								}

								if((ext-extTprev)==200){
										if(asicNo0==asicNo){
												addflag=true;
										}else if(!nextflag){
												nextflag=true; 
												inext=i;
										}
								}else if((ext-extTprev)>200){ 
										endflag=true;
								}else if((ext-extTprev)>=0 && (ext-extTprev)<200){
										if(!nextflag){
												nextflag=true; 
												inext=i;
										}
								}else{
										if(!rewind_flag || (r_i%1000==0)){
												cout<<"The re-written event seemed started "<<ext<<" and \t"<<extTprev<<" is delta extT and extTprev at "
														<<iprev<<endl<<i<<"'th entry\t"<<r_i<<" entries rewinded."<<" The remaining entries are "<<n-i<<endl;

												rewind_flag=true;
										}
										r_i++;
								}

						}

				}

				if(extT0 < extT0prev && (n-i)<10000 ) {
						break;
				}else if(extT0 <extT0prev && (n-i)>10000) {
						cout<<"It seems not just the rewinded entries"<<endl;
						return -1;
				}

				if(i0%100000==0) { std::cout<<"\r"<<i<<"entry done. "<<"timespent : "<<time(&tm)-time0<<"s. The tmap2 size is "<<tmap2.size(); std::cout.flush();}
				if(nextflag) i=inext;
				i0++;

		}
		std::cout<<endl<<"ASIC multiplexer timestamp shift mapping done.\n Calibration applying and filling tree start"<<endl;

		//		time0=time(&tm);
		//		itmap=tmap.begin();
		//		itmap2=tmap2.begin();	
		//		it_indexer=indexer.begin();

		if((tmap2.size() <= n) && (tmap.size()==(n-starti))){ cout<<"The sizes of the all C-structures are fine"<<endl;}
		else {cout<<"The sizes are not good"<<tmap.size()<<"\t"<<tmap2.size()<<"\t"<<"\t"<<n<<endl; return -1;}

		for(i=starti;i<tmap2.size();i++){

				if(i%100000==0) {std::cout<<"\r"<<i<<"'s entry done. Timespent "<<time(&tm)-time0<<"s";std::cout.flush();}

				//				j=indexer[i];
				j=tmap2[i-starti].first;
				if(j<0) continue;

				tr->GetEntry(j);
				old_ext=tmap[j-starti];
				//				itmap2=tmap2.find(j);
				//				ext=itmap2->second;
				ext=tmap2[i-starti].second;
				//				it_indexer++;

				ct->GetEntry((dssdNo)*256+stripNo);

				if(quality==3) continue;

				if(infoCode==0 && feeNo<25 && feeNo>0){
						if(rangeType==0) {
								//R59 setting
								//if((adcData>thresh && thresh>-500) || thresh<-400)
								//								if(adcData>thresh)
								//								{
								//										if(quality>1) {
								/*					if(quality==4){
													adcData_off=(-offset+adcData);					//with gain correction
													}

													adcData_off=(-offset+adcData)*0.158/gain;					//with gain correction
								 */

								adcData_off=(-offset+adcData);					//without gain correction
								//												quality=1;
								//												newt->Fill();
								//										}

								/*	if(stripNo<128){
									/											adcData_off= (-offset + adcData)*0.013/gain;
									}else if(stripNo>=128){
									adcData_off= (-offset + adcData)*0.077/gain;
									}												//R18 setting
									if(adcData_off < thresh) quality=-1;
									if(quality==1)	newt->Fill();
								 */


								//				if(quality==1 && adcData_off >1500) quality=0;		//check for noise petestal suspect (20170103)
								//								}
						}
						else if(rangeType==1) {
								//								quality=1;								// For auto calibration, adcData of any quality are accepted
								adcData_off= adcData;
								//								newt->Fill();
						}

						newt->Fill();

				}

		}
		std::cout<<std::endl;
		nf->WriteTObject(newt);
		nf->Close();

		return 1;
}
