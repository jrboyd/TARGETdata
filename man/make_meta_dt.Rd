% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/append.R
\name{make_meta_dt}
\alias{make_meta_dt}
\alias{make_meta_dt.full}
\title{make_meta_dt}
\usage{
make_meta_dt(
  dt,
  clin_dt = load_clinical_data(),
  splice_clust_dt = load_splice_clusters(),
  include_all = FALSE,
  extra_vars_included = character()
)

make_meta_dt.full(dt, clin_dt = load_clinical_data())
}
\arguments{
\item{dt}{data.table or matrix from load_*.}

\item{clin_dt}{clinical info, as from load_clinical_data().}

\item{include_all}{if TRUE all clinical attributes are included.}
}
\value{
data.table of clinical data fro all sample_ids
}
\description{
make_meta_dt

make_meta_dt.full
}
\examples{
mat = load_miR_RPM()
mat = filter_expression_to_valid(mat)
meta_dt = make_meta_dt(mat)

dt = convert_matrix_2_data.table(mat)
make_meta_dt(dt)

make_meta_dt.full(mat)
make_meta_dt.full(dt)

}
