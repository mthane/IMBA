#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>

template<typename T>
class partition_generator
{
  public:
    size_t m;
    vector<vector<vector<T> > > RES;
    vector<T> &L;
    vector<int> A;
    size_t s;
    partition_generator(std::vector<T> &LD,size_t mb):
      m(mb),
      L(LD),
      s(0)
  {
    A=vector<int>(L.size());
  }
    void generate_subsets()
    {
      vector<vector<T> > ss=vector<vector<T> > (m);
      for(auto i=0;i<L.size();i++)
      {
        ss[A[i]].push_back(L[i]);
      }
      //std::cout << "A: " << printVector<int>(A) << std::endl;
      RES.push_back(ss);
      s++;
    }

    void m_partitions_of_n()
    {
      size_t n=L.size();
      for(auto j=0;j<m;++j)
        A[n-m+j]=j;

      f(m,n,0);
      std::vector<std::vector<std::vector<T> > > ADD;
      for(auto &i:RES)
      {
        std::vector<int> pbase;
        for(auto k=0;k<m;k++)
          pbase.push_back(k);

        std::next_permutation(pbase.begin(),pbase.end());

        do {
          vector<vector<T> > ss=vector<vector<T> > (m);
          for(auto k=0;k<m;k++)
            ss[k]=i[pbase[k]];

          ADD.push_back(ss);
        } while ( std::next_permutation(pbase.begin(),pbase.end()) );
      }

     size_t sz = RES.size();

     RES.insert(RES.end(),
         std::make_move_iterator(ADD.begin()),
         std::make_move_iterator(ADD.end()));

     std::inplace_merge(RES.begin(), RES.begin() + sz, RES.end()); 
    }

    void f(size_t mu, size_t nu, size_t sigma)
    {
      if(mu==2)
        generate_subsets();
      else
        f(mu-1,nu-1,(mu+sigma)%2);

      if(nu==mu+1)
      {
        A[mu-1]=mu-1;
        generate_subsets();
        while(A[nu-1]>0)
        {
          A[nu-1]=A[nu-1]-1;
          generate_subsets();
        }
      }
      else if(nu>mu+1)
      {
        if((mu+sigma)%2==1)
          A[nu-2]=mu-1;
        else
          A[mu-1]=mu-1;

        if((A[nu-1]+sigma)%2==1)
          b(mu,nu-1,0);
        else
          f(mu,nu-1,0);
          
        while(A[nu-1]>0)
        {
          A[nu-1]=A[nu-1]-1;
          if((A[nu-1]+sigma)%2==1)
            b(mu,nu-1,0);
          else
            f(mu,nu-1,0);
        }
      }
    }

    void b(size_t mu, size_t nu, size_t sigma)
    {
      if(nu==mu+1)
      {
        while(A[nu-1]<mu-1)
        {
          generate_subsets();
          ++A[nu-1];
        }
        generate_subsets();
        A[mu-1]=0;
      }
      else if(nu>mu+1)
      {
        if((A[nu-1]+sigma)%2==1)
          f(mu,nu-1,0);
        else
          b(mu,nu-1,0);
        while(A[nu-1]<mu-1)
        {
          ++A[nu-1];
          if((A[nu-1]+sigma)%2==1)
            f(mu,nu-1,0);
          else
            b(mu,nu-1,0);
        }
        if((mu+sigma)%2==1)
          A[nu-2]=0;
        else
          A[mu-1]=0;
      }
      if(mu==2)
        generate_subsets();
      else
        b(mu-1,nu-1,(mu+sigma)%2);
    }

};
