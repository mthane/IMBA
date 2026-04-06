#include <opencv2/core.hpp>
#include <lpsolve/lp_lib.h>
#undef REAL
#undef NORMAL
#include <string>
#include <iostream>
#include <vector>
#include <stack>
#include <map>
#include "larvaObject.hpp"
#include "lrvTrackDebug.hpp"
//#include <SWI-Prolog.h>
//#include <SWI-cpp.h>

using namespace cv;
using namespace std;
using namespace lrvTrack;

extern std::map<size_t, std::vector<size_t> > parent_blobs;
extern std::map<size_t, std::vector<size_t> > children_blobs;
extern std::map<size_t,larvaObject> detected_larvae;
extern double LRVTRACK_MPP;
/*bool testPLQuery()
{
  char *expression=strdup("library(clpfd)");
  char *plav[2];
  char *str1=strdup("-g use_module(library(clpfd)).");
  //char *str1=strdup("");
  plav[0] = str1;
  plav[1] = NULL;
  if ( !PL_initialise(1, plav) )
  {
        PL_halt(1);
        return false;
  }
  else
  {
    PlTermv av(1);
    PlTermv ava("library(clpfd)");

    PlQuery w("use_module", ava);
    while( w.next_solution() )
    {
    }
    PlQuery q("current_module", av);
    while( q.next_solution() )
      cout << (char *)av[0] << endl;

        return TRUE;

  }
  return true;
}*/

bool TUmatrixTest(Mat m)
{
  //cout << "Whole matrix: " << endl << m << endl;
  size_t minDim=min(m.rows,m.cols);
  for(size_t i=2;i<=minDim;i++)
  {
    for(size_t posx=0;posx<=(m.cols-i);posx++)
    {
      for(size_t posy=0;posy<=(m.rows-i);posy++)
      {
        Mat s=m(Range(posy,posy+i),Range(posx,posx+i));
        //cout << "Testing matrix: " << endl << s << endl;
        if(fabs(determinant(s))>1)
        {
          return false;
        }
      }
    }
  }
  return true;
}

void LPmaxFunction(Mat &mf)
{
  vector<double> mf_vec;
  //TODO 10000 pixels should be enough but we should perhaps standardize
  for (auto &p: detected_larvae)
  {
    larvaObject &l=p.second;
    mf_vec.push_back(-(1.0-(l.area_mean/10000.0)));
  }
  Mat m(mf_vec);
  m.copyTo(mf);
}
int LPconstrAvg(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex,
	    double avgSize)
{
  size_t total_larvae=detected_larvae.size();
  int ret;
  //The vector containing the contraints to be converted to Mat
  //List of nodes that have been scanned as parents either directly
  //or through the deep search loop.
  vector<vector<double> > constraints;
  vector<size_t> parents_scanned;
  /*//List of nodes that converge to themselves in two steps:
  [>           
   *  \ /        A
       |        / \
       |       B   C
       |        \ /
       V         D
  <]
  vector<size_t> marked_nodes;
  size_t idx=0;
  for(auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;

    indexToID.push_back(lID);
    IDtoIndex[lID]=idx;
    idx++;
    //If all children of current node have only one parent
    // HERE we assume that these assignments must be correct
    // and if a parent has a child then that child has him as parent
    if((all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
            [](size_t i){return 
            (parent_blobs[i].size()==1);})) && 
        //If there's more than one children
        children_blobs[lID].size()>1 &&
        //If all children have only one parent
        (all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
                [](size_t i){return 
                (children_blobs[i].size()==1);})) )
    {
      //Check if they all converge to the same node
      size_t ccID=children_blobs[children_blobs[lID].front()].front();
      class ccomp{
        public:
          ccomp(size_t a):
            comp_cID(a){}
          bool converges(size_t a){return children_blobs[a].front()==comp_cID;}
        private:
          size_t comp_cID;
      };
      ccomp cur_ccomp(ccID);
      //See here for explanation:
      //http://stackoverflow.com/a/6854152
      if(all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
            std::bind(&ccomp::converges, &cur_ccomp,std::placeholders::_1)
            )
        )
      {
        marked_nodes.push_back(lID);
      }
    }
  }*/

  /*if(testPLQuery())
    cout << "PL TRUE" << endl;*/
  lprec *lp;
  lp = make_lp(0, detected_larvae.size());
  set_add_rowmode(lp, TRUE);
  size_t vars=detected_larvae.size();
  vector<int> cols;
  for(size_t i=0;i<vars;i++)
    cols.push_back(i+1);
  int  *colno=&cols[0];
  double minarea=DBL_MAX;
  for(auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    //cout << "CONSTRAINT_CONSTRUCT: P[" << lID << "]=" << printVector(parent_blobs[lID]) << endl;
    //cout << "CONSTRAINT_CONSTRUCT: C[" << lID << "]=" << printVector(children_blobs[lID]) << endl;
    size_t minVal=max(parent_blobs[lID].size(),
        children_blobs[lID].size());
    if(minVal==0)
      minVal=1;
    int con_type=GE;
    if(p.second.area_mean>=2*avgSize)
    {
      con_type=GE;
      minVal=2;
    }
    if(p.second.area_mean<=0.9*avgSize)
    {
      con_type=EQ;
      minVal=1;
    }

    vector<double> minConst(total_larvae,0.0);
    minConst[IDtoIndex[lID]]=1.0;
    //minConst.push_back(minVal);
    double* row = &minConst[0];
    if(!add_constraintex(lp, vars, row, colno, con_type, minVal))
      return(3);
    if(detected_larvae[lID].area_mean<minarea)
      minarea=detected_larvae[lID].area_mean<minarea;
    /*for(auto &c: children_blobs[lID])
    {
      for(auto &d: children_blobs[lID])
      {
        if(c!=d)
        {
          if(detected_larvae[c].area_mean>
             detected_larvae[d].area_mean)
          {
            vector<double> leConst(total_larvae,0.0);
            leConst[IDtoIndex[c]]=1.0;
            leConst[IDtoIndex[d]]=-1.0;
            double* lerow = &leConst[0];
            if(!add_constraintex(lp, vars, lerow, colno, GE, 0.0))
              return(3);
          }
        }
      }
    }*/
  }
  /*for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    vector<double> constr(total_larvae,0.0);
    if(detected_larvae[lID].area_mean<minarea*1.5)
    {
      constr[IDtoIndex[lID]]=1.0;
      double* row = &constr[0];
      if(!add_constraintex(lp, vars, row, colno, EQ,1.0))
        return(3);
    }
  }*/
  int col=0;
  for (auto &p: detected_larvae)
  {
    set_int(lp,col++,TRUE);
    size_t lID=p.second.larva_ID;
    //If this parent is in the parents_scanned, skip
    if(find(parents_scanned.begin(),parents_scanned.end(),lID)
        !=parents_scanned.end())
      continue;

    vector<double> constr(total_larvae,0.0);
    //vector<double> rconstr(total_larvae,0.0);

    //Stack of parents for the search of parents sharing children
    //with the current node.
    stack<size_t> sParents;
    vector<size_t> cChildren;
    vector<size_t> cParents;
    sParents.push(lID);
    while(!sParents.empty())
    {
      size_t pID=sParents.top();
      sParents.pop();
      parents_scanned.push_back(pID);
      cParents.push_back(pID);
      
      for(auto &child: children_blobs[pID])
      {
        if(find(cChildren.begin(),cChildren.end(),child)
            ==cChildren.end())
        {
          cChildren.push_back(child);
          for(auto &parent: parent_blobs[child])
          {
            if(find(cParents.begin(),cParents.end(),parent)
                ==cParents.end())
            {
              sParents.push(parent);
            }
          }
        }
      }
    }
    for(auto &p:cParents)
    {
      constr[IDtoIndex[p]]=-1.0;
      //rconstr[IDtoIndex[p]]=1.0;
    }
    for(auto &c:cChildren)
    {
      constr[IDtoIndex[c]]=1.0;
      //rconstr[IDtoIndex[c]]=-1.0;
    }
    //constr.push_back(0);
    //rconstr.push_back(0);


    if(children_blobs[lID].size()!=0)
    {
      double* row = &constr[0];
      if(!add_constraintex(lp, vars, row, colno, EQ, 0.0))
        return(3);
    }
  }
  
  set_add_rowmode(lp, FALSE);

  vector<double> mf_vec;
  for (auto &p: detected_larvae)
  {
    larvaObject &l=p.second;
    mf_vec.push_back(1.0-(l.area_mean/20000.0)); // 20000 -> MAX_OBJ_SIZE as defined in header??
  }
  double* row = &mf_vec[0];
  /* set the objective in lpsolve */
  if(!set_obj_fnex(lp, vars, row, colno))
    return(4);

  set_minim(lp);
  write_LP(lp, stdout);
  /* write_lp(lp, "model.lp"); */

  /* I only want to see important messages on screen while solving */
  set_verbose(lp, IMPORTANT);
  ret = solve(lp);
  if(ret == OPTIMAL)
  {
    ret = 0;
    cout << "LPSOLVE: Found " << get_solutioncount(lp) << " solutions." << endl;
  }
  else
  {
    cerr << "LPSOLVE: Couldn't solve system" << ret << endl;
    return(5);
  }
  /* a solution was calculated, now lets get some results */
  /* variable values */

  double *results = (double *) malloc(vars * sizeof(double));
  get_variables(lp, results);
  for(size_t i=0;i<vars;i++)
    solution.push_back(static_cast<int>(results[i]+0.5));

  cerr << "LP Solution : [";
  for(auto &v:solution)
    cerr << v << " ";
  cerr << "]" <<endl;
  /* free allocated memory */
  if(results != NULL)
    free(results);
  if(lp != NULL) {
    /* clean up such that all used memory by lpsolve is freed */
    delete_lp(lp);
  }
  return 0;

  /*Mat MatConstr=Mat(constraints.size(),total_larvae+1,CV_64FC1);
  for(int r=0;r<MatConstr.rows;r++)
  {
    Mat crow(constraints[r]);
    crow=crow.t();
    crow.copyTo(MatConstr.row(r));
  }
  MatConstr.copyTo(mConstr);
  [>if(!TUmatrixTest(mConstr(Range(0,mConstr.rows),Range(0,mConstr.cols-1))))
    BOOST_LOG_TRIVIAL(debug) << "CONSTRAINTS: WARN!! : Constraints matrix is not TU" << endl;
  else
    BOOST_LOG_TRIVIAL(debug) << "CONSTRAINTS: Constraints matrix is TU" << endl;<]*/
}

//Construct constraints for Simplex:
//Ax<=b
//A constraints of the form        
//-------------- A ---------------  b
//0  -1  1   -1 0  0  0  0  0  ...  0.1
//p1 p2  p3  p4 p5 p6 p7 p8 p9 ... 
//0p1 -p2 + p3 -p4  + 0p5 ... 
//and the oposite:
//0  1   -1   1  0  0  0  0  0 ... 0
//p1 p2  p3  p4 p5 p6 p7 p8 p9 ...
//0p1 +p2 -p3 +p4  + 0p5 ...
int LPconstr(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex)
{
  size_t total_larvae=detected_larvae.size();
  int ret;
  //The vector containing the contraints to be converted to Mat
  //List of nodes that have been scanned as parents either directly
  //or through the deep search loop.
  vector<vector<double> > constraints;
  vector<size_t> parents_scanned;
  /*//List of nodes that converge to themselves in two steps:
  [>           
   *  \ /        A
       |        / \
       |       B   C
       |        \ /
       V         D
  <]
  vector<size_t> marked_nodes;
  size_t idx=0;
  for(auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;

    indexToID.push_back(lID);
    IDtoIndex[lID]=idx;
    idx++;
    //If all children of current node have only one parent
    // HERE we assume that these assignments must be correct
    // and if a parent has a child then that child has him as parent
    if((all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
            [](size_t i){return 
            (parent_blobs[i].size()==1);})) && 
        //If there's more than one children
        children_blobs[lID].size()>1 &&
        //If all children have only one parent
        (all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
                [](size_t i){return 
                (children_blobs[i].size()==1);})) )
    {
      //Check if they all converge to the same node
      size_t ccID=children_blobs[children_blobs[lID].front()].front();
      class ccomp{
        public:
          ccomp(size_t a):
            comp_cID(a){}
          bool converges(size_t a){return children_blobs[a].front()==comp_cID;}
        private:
          size_t comp_cID;
      };
      ccomp cur_ccomp(ccID);
      //See here for explanation:
      //http://stackoverflow.com/a/6854152
      if(all_of(children_blobs[lID].begin(),children_blobs[lID].end(),
            std::bind(&ccomp::converges, &cur_ccomp,std::placeholders::_1)
            )
        )
      {
        marked_nodes.push_back(lID);
      }
    }
  }*/

  /*if(testPLQuery())
    cout << "PL TRUE" << endl;*/
  lprec *lp;
  lp = make_lp(0, detected_larvae.size());
  set_add_rowmode(lp, TRUE);
  size_t vars=detected_larvae.size();
  vector<int> cols;
  for(size_t i=0;i<vars;i++)
    cols.push_back(i+1);
  int  *colno=&cols[0];
  double minarea=DBL_MAX;
  for(auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    //cout << "CONSTRAINT_CONSTRUCT: P[" << lID << "]=" << printVector(parent_blobs[lID]) << endl;
    //cout << "CONSTRAINT_CONSTRUCT: C[" << lID << "]=" << printVector(children_blobs[lID]) << endl;
    size_t minVal=max(parent_blobs[lID].size(),
        children_blobs[lID].size());
    if(minVal==0)
      minVal=1;
    int con_type=GE; // Ax >= b
    /*if(p.second.area_mean*LRVTRACK_MPP*LRVTRACK_MPP>=3.95)
    {
      con_type=GE;
      minVal=2;
    }*/
    /*if(p.second.area_mean*LRVTRACK_MPP*LRVTRACK_MPP<=2.95)
    {
      con_type=EQ;
      minVal=1;
    }*/

    vector<double> minConst(total_larvae,0.0);
    minConst[IDtoIndex[lID]]=1.0; // jede ID hat maximal 1 index
    //minConst.push_back(minVal);
    double* row = &minConst[0];
// Ax >= b?
// 0*x1 + 1*x2 + 0*x3+... + 0*xn >= minVal
    if(!add_constraintex(lp, vars, row, colno, con_type, minVal))
      return(3);
    if(detected_larvae[lID].area_mean<minarea)
      minarea=detected_larvae[lID].area_mean<minarea;
    /*for(auto &c: children_blobs[lID])
    {
      for(auto &d: children_blobs[lID])
      {
        if(c!=d)
        {
          if(detected_larvae[c].area_mean>
             detected_larvae[d].area_mean)
          {
            vector<double> leConst(total_larvae,0.0);
            leConst[IDtoIndex[c]]=1.0;
            leConst[IDtoIndex[d]]=-1.0;
            double* lerow = &leConst[0];
            if(!add_constraintex(lp, vars, lerow, colno, GE, 0.0))
              return(3);
          }
        }
      }
    }*/
  }
  /*for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    vector<double> constr(total_larvae,0.0);
    if(detected_larvae[lID].area_mean<minarea*1.5)
    {
      constr[IDtoIndex[lID]]=1.0;
      double* row = &constr[0];
      if(!add_constraintex(lp, vars, row, colno, EQ,1.0))
        return(3);
    }
  }*/
  int col=0;
  for (auto &p: detected_larvae)
  {
    set_int(lp,col++,TRUE);
    size_t lID=p.second.larva_ID;
    //If this parent is in the parents_scanned, skip
    if(find(parents_scanned.begin(),parents_scanned.end(),lID)
        !=parents_scanned.end())
      continue;

    vector<double> constr(total_larvae,0.0);
    //vector<double> rconstr(total_larvae,0.0);

    //Stack of parents for the search of parents sharing children
    //with the current node.
    stack<size_t> sParents;
    vector<size_t> cChildren;
    vector<size_t> cParents;
    sParents.push(lID);
    while(!sParents.empty())
    {
      size_t pID=sParents.top();
      sParents.pop();
      parents_scanned.push_back(pID);
      cParents.push_back(pID);
      
      for(auto &child: children_blobs[pID])
      {
        if(find(cChildren.begin(),cChildren.end(),child)
            ==cChildren.end())
        {
          cChildren.push_back(child);
          for(auto &parent: parent_blobs[child])
          {
            if(find(cParents.begin(),cParents.end(),parent)
                ==cParents.end())
            {
              sParents.push(parent);
            }
          }
        }
      }
    }
    for(auto &p:cParents)
    {
      constr[IDtoIndex[p]]=-1.0;
      //rconstr[IDtoIndex[p]]=1.0;
    }
    for(auto &c:cChildren)
    {
      constr[IDtoIndex[c]]=1.0;
      //rconstr[IDtoIndex[c]]=-1.0;
    }
    //constr.push_back(0);
    //rconstr.push_back(0);


    if(children_blobs[lID].size()!=0)
    {
      double* row = &constr[0];
      if(!add_constraintex(lp, vars, row, colno, EQ, 0.0))
        return(3);
    }
  }
  
  set_add_rowmode(lp, FALSE);

  vector<double> mf_vec;
  for (auto &p: detected_larvae)
  {
    larvaObject &l=p.second;
    mf_vec.push_back(1.0-(l.area_mean/20000.0));
  }
  double* row = &mf_vec[0];
  /* set the objective in lpsolve */
  if(!set_obj_fnex(lp, vars, row, colno))
    return(4);

  set_minim(lp);
  write_LP(lp, stdout);
  /* write_lp(lp, "model.lp"); */

  /* I only want to see important messages on screen while solving */
  set_verbose(lp, IMPORTANT);
  ret = solve(lp);
  if(ret == OPTIMAL)
  {
    ret = 0;
    cout << "LPSOLVE: Found " << get_solutioncount(lp) << " solutions." << endl;
  }
  else
  {
    cerr << "LPSOLVE: Couldn't solve system" << ret << endl;
    return(5);
  }
  /* a solution was calculated, now lets get some results */
  /* variable values */

  double *results = (double *) malloc(vars * sizeof(double));
  get_variables(lp, results);
  for(size_t i=0;i<vars;i++)
    solution.push_back(static_cast<int>(results[i]+0.5));

  cerr << "LP Solution : [";
  for(auto &v:solution)
    cerr << v << " ";
  cerr << "]" <<endl;
  /* free allocated memory */
  if(results != NULL)
    free(results);
  if(lp != NULL) {
    /* clean up such that all used memory by lpsolve is freed */
    delete_lp(lp);
  }
  return 0;

  /*Mat MatConstr=Mat(constraints.size(),total_larvae+1,CV_64FC1);
  for(int r=0;r<MatConstr.rows;r++)
  {
    Mat crow(constraints[r]);
    crow=crow.t();
    crow.copyTo(MatConstr.row(r));
  }
  MatConstr.copyTo(mConstr);
  [>if(!TUmatrixTest(mConstr(Range(0,mConstr.rows),Range(0,mConstr.cols-1))))
    BOOST_LOG_TRIVIAL(debug) << "CONSTRAINTS: WARN!! : Constraints matrix is not TU" << endl;
  else
    BOOST_LOG_TRIVIAL(debug) << "CONSTRAINTS: Constraints matrix is TU" << endl;<]*/
}

//Construct constraints for Simplex:
//Ax<=b
//A constraints of the form        
//-------------- A ---------------  b
//0  -1  1   -1 0  0  0  0  0  ...  0.1
//p1 p2  p3  p4 p5 p6 p7 p8 p9 ... 
//0p1 -p2 + p3 -p4  + 0p5 ... 
//and the oposite:
//0  1   -1   1  0  0  0  0  0 ... 0
//p1 p2  p3  p4 p5 p6 p7 p8 p9 ...
//0p1 +p2 -p3 +p4  + 0p5 ...
int LPconstr(std::vector<double> &solution,
             std::vector<size_t> &indexToID, 
            std::map<size_t,size_t> &IDtoIndex,
            double avg_size)
{
  size_t total_larvae=detected_larvae.size();
  int ret;
  //The vector containing the contraints to be converted to Mat
  //List of nodes that have been scanned as parents either directly
  //or through the deep search loop.
  vector<vector<double> > constraints;

  lprec *lp;
  lp = make_lp(0, total_larvae*2);
  vector<int> cols;
  for(size_t i=0;i<2*total_larvae;i++)
    cols.push_back(i+1);
  int  *colno=&cols[0];

  vector<double> mf_vec;
  for (auto &p: detected_larvae)
  {
    mf_vec.push_back(0.0);
  }
  for (auto &p: detected_larvae)
  {
    mf_vec.push_back(1.0);
  }
  double* row = &mf_vec[0];
  /* set the objective in lpsolve */
  if(!set_obj_fnex(lp, 2*total_larvae, row, colno))
    return(4);

  set_add_rowmode(lp, TRUE);
  //Constraints of minimum larvae in blob
  double minarea=DBL_MAX;
  for(auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    if(detected_larvae[lID].area_mean<minarea)
      minarea=detected_larvae[lID].area_mean<minarea;
    //cout << "CONSTRAINT_CONSTRUCT: P[" << lID << "]=" 
    //  << printVector(parent_blobs[lID]) << endl;
    //cout << "CONSTRAINT_CONSTRUCT: C[" << lID << "]=" 
    //  << printVector(children_blobs[lID]) << endl;
    size_t minVal=max(parent_blobs[lID].size(),
        children_blobs[lID].size());
    if(minVal==0)
      minVal=1;
    vector<double> minConst(2*total_larvae,0.0);
    minConst[IDtoIndex[lID]]=1.0;
    //minConst.push_back(minVal);
    double* row = &minConst[0];
    if(!add_constraintex(lp, 2*total_larvae, row, colno, GE, minVal))
      return(3);
    for(auto &c: children_blobs[lID])
    {
      for(auto &d: children_blobs[lID])
      {
        if(c!=d)
        {
          if(detected_larvae[c].area_mean>
             detected_larvae[d].area_mean)
          {
            vector<double> leConst(total_larvae,0.0);
            leConst[IDtoIndex[c]]=1.0;
            leConst[IDtoIndex[d]]=-1.0;
            double* lerow = &leConst[0];
            if(!add_constraintex(lp, total_larvae, lerow, colno, GE, 0.0))
              return(3);
          }
        }
      }
    }
  }

  for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    vector<double> constr(2*total_larvae,0.0);
    if(detected_larvae[lID].area_mean<minarea*1.5)
    {
      constr[IDtoIndex[lID]]=1.0;
      double* row = &constr[0];
      if(!add_constraintex(lp, 2*total_larvae, row, colno, EQ,1.0))
        return(3);
    }
  }

  vector<size_t> parents_scanned;
  for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    //If this parent is in the parents_scanned, skip
    if(find(parents_scanned.begin(),parents_scanned.end(),lID)
        !=parents_scanned.end())
      continue;

    vector<double> constr(2*total_larvae+1,0.0);

    //Stack of parents for the search of parents sharing children
    //with the current node.
    stack<size_t> sParents;
    vector<size_t> cChildren;
    vector<size_t> cParents;
    sParents.push(lID);
    while(!sParents.empty())
    {
      size_t pID=sParents.top();
      sParents.pop();
      parents_scanned.push_back(pID);
      cParents.push_back(pID);
      
      for(auto &child: children_blobs[pID])
      {
        if(find(cChildren.begin(),cChildren.end(),child)
            ==cChildren.end())
        {
          cChildren.push_back(child);
          for(auto &parent: parent_blobs[child])
          {
            if(find(cParents.begin(),cParents.end(),parent)
                ==cParents.end())
            {
              sParents.push(parent);
            }
          }
        }
      }
    }
    for(auto &p:cParents)
    {
      constr[IDtoIndex[p]]=-1.0;
      //rconstr[IDtoIndex[p]]=1.0;
    }
    for(auto &c:cChildren)
    {
      constr[IDtoIndex[c]]=1.0;
      //rconstr[IDtoIndex[c]]=-1.0;
    }
    //constr.push_back(0);
    //rconstr.push_back(0);


    if(children_blobs[lID].size()!=0)
    {
      double* row = &constr[0];
      if(!add_constraintex(lp, 2*total_larvae, row, colno, EQ, 0.0))
        return(3);
    }
  }

  //Construct the absolute constraints...
  //We want our objective funtion to minimize:
  //           | Size_1 - Count_1*avg_size | + ... + | Size_n - Count_n*avg_size |
  //           Size is the actual area of the blob and count_n 
  //           are the contents we are trying to find.
  //           We need two constraints:
  //            X'_i>= Size_i - Count_i * avg_size,
  //            X'_i>= Count_i * avg_size - Size_i
  //           Then the objective function is min: X'_1+...+X'_n
  for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    //double d=(p.second.end_frame-p.second.start_frame)/1000;
    double d=1.0;
    vector<double> absConst1(2*total_larvae,0.0);
    vector<double> absConst2(2*total_larvae,0.0);
    absConst1[IDtoIndex[lID]]=d*avg_size;
    absConst1[IDtoIndex[lID]+total_larvae]=1.0;
    absConst2[IDtoIndex[lID]]=-d*avg_size;
    absConst2[IDtoIndex[lID]+total_larvae]=1.0;
    double* abrow1 = &absConst1[0];
    double* abrow2 = &absConst2[0];
      
    if(!add_constraintex(lp, 2*total_larvae, abrow1, colno, GE, d*p.second.area_mean))
      return(3);
    if(!add_constraintex(lp, 2*total_larvae, abrow2, colno, GE, -d*p.second.area_mean))
      return(3);
  
  }
  set_add_rowmode(lp, FALSE);

  set_minim(lp);
  write_LP(lp, stdout);
  /* write_lp(lp, "model.lp"); */

  /* I only want to see important messages on screen while solving */
  set_verbose(lp, DETAILED);
  ret = solve(lp);
  if(ret == OPTIMAL)
    ret = 0;
  else
  {
    cerr << "LPSOLVE: Couldn't solve system: " << ret << endl;
    return(5);
  }

  /* a solution was calculated, now lets get some results */
  /* variable values */

  double *results = (double *) malloc(2*total_larvae* sizeof(double));
  get_variables(lp, results);
  for(size_t i=0;i<2*total_larvae;i++)
    solution.push_back(static_cast<int>(results[i]+0.5));

  cerr << "LP Solution : [";
  for(auto &v:solution)
    cerr << v << " ";
  cerr << "]" <<endl;
  /* free allocated memory */
  if(results != NULL)
    free(results);
  if(lp != NULL) {
    /* clean up such that all used memory by lpsolve is freed */
    delete_lp(lp);
  }
  return 0;
}

