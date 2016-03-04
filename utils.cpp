//RG:WIP NO NEED FOR lower_bound AND upper_bound. WHOLE ARRAY IS PASSED
int bisection_search(doub target, doub arr[], doub lower_bound, doub upper_bound, int size, int& target_idx, int& lower_bound_idx, int& upper_bound_idx) { 
  // Given array 'arr' and target value 'target' return index 'idx' where arr[idx]=target

  doub closest_to_target;
  // doub below_target=arr[lower_bound];
  // doub above_target=arr[upper_bound];
  doub below_target=lower_bound;
  doub above_target=upper_bound;

  if((target<=above_target) && (below_target<=target)){

	while(upper_bound_idx>lower_bound_idx+1){ // [lower_bound_idx,upper_bound_idx]
      target_idx=(lower_bound_idx+upper_bound_idx)>>1; // RG: bit-shift operations for speed up? (lower_bound_idx+upper_bound_idx)>>1 = (lower_bound_idx+upper_bound_idx)/2 ? golden ratio is faster than 0.5 interval (optimal 1st order method), even better Brent's method
      closest_to_target=arr[target_idx];
      if(target<closest_to_target){
        upper_bound_idx=target_idx; // [lower_bound_idx,target_idx]
      } else {
        lower_bound_idx=target_idx; // [target_idx,upper_bound_idx]
      };
	};

  } else {
    printf(YELLOW"[evalpointzero.cpp]:"RED" Lookup error of closest geodesic point below_target=%f target=%f above_target=%f\nEXITING\n"RESET,below_target,target,above_target);
  //t_solvetrans += (clock() - t_b4_solvetrans) / (doub)CLOCKS_PER_SEC;
    exit(-1);
  };

  return (0);
}



/******************************/
/* RG:OUTPUT TEMPERATURE INFO */

void temperature_diag(doub rr, doub costh, doub ph, doub T_sim, doub tet, doub tpt, doub rho, doub rhonor, int currth, int geodesic_idx) { 
  // "OUTPUT TEMPERATURE INFO TO FILE FOR DIAGNOSTIC PURPOSES"
  // ALONG GEODESIC

  // printf(YELLOW"[evalpointzero.cpp]: "RESET"TEMPERATURE DIAGNOSTIC!\n");

  stringstream temperature_diag_sstr;
  FILE * TeTp_file; 

  temperature_diag_sstr<<"temperature_diag"<<(int)100*a<<"th"<<(int)floor(100*th+1e-6)<<"fn"<<fnum;
  string append_label;
  append_label=temperature_diag_sstr.str();

  // if (iiy%geodesic_output_every_x==0 && iix==nxy/2) { // sample x=const slice of picture plane
  // if (true) { // unlimited output...
  if (false) {
  // if (choose_geodesic==500) { // limited output // RG:ERROR geodesic_idx is the number of points on a geodesic: Why do I get output at all?!
  // if (curr==1) { // limited output
    
    // int geodesic_label=ix;
    stringstream temperature_diag_sstr_append;
    temperature_diag_sstr_append<<append_label<<"geod"<<geodesic_idx;
    string stra = temperature_diag_sstr_append.str();
    
    //RG: should really open/close file outside computation...
    TeTp_file=fopen ((dir+stra+".dat").c_str(),"a");
    
    //for (int geo_idx=0; geo_idx<=maxco; geo_idx++) { // maxco
    // for (int geo_idx=0; geo_idx<=stN; geo_idx++) { // maxco
      // for (int p=0; p<=7; p++) { // r:p-> costh: ph:
      //   // fprintf(TeTp_file,"%f ",ppy[currth].cooxx[p][geo_idx]);
      //   fprintf(TeTp_file,"%f ",ppy[curr].cooxx[p][geodesic_idx]);
      // }
      // lamx[] vs lam[]?
      // fprintf(TeTp_file,"%f %d %d\n",ppy[currth].lamx[geo_idx],geodesic_label,geo_idx);
    static int one_time_only = 0;
    if (one_time_only==0) {
      fprintf(TeTp_file,"# rr \t costh \t   ph \t    tsim     t_e_mod  tpt      t_e_orig     rho\n");
      one_time_only++;
    }
    fprintf(TeTp_file,"%.2e %.2e %.2e %.2e %.2e %.2e %.2e\n",rr,costh,ph,T_sim,tet,tpt,rho/rhonor); // ts[ia] ~ ts[ib]
    
    fclose(TeTp_file);
    }

    // } //for (int geodesic_label=0; geo_idx<=maxco; geo_idx+=1000) { // if (curr)

} // if TEMPERATURE_DIAGNOSTIC
