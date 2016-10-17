struct NRsparseMat
{
	Int nrows;
	Int ncols;
	Int nvals;
	VecInt col_ptr;
	VecInt row_ind;
	VecDoub val;

	NRsparseMat(Int m,Int n,Int nnvals);
};
NRsparseMat::NRsparseMat(Int m,Int n,Int nnvals) : nrows(m),ncols(n),
	nvals(nnvals),col_ptr(n+1,0),row_ind(nnvals,0),val(nnvals,0.0) {}

