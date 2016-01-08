# echo "cleanup GeoPrimitives GeoPrimitives-00-00-33-05 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/cern.ch/sw/contrib/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtGeoPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtGeoPrimitivestempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtGeoPrimitivestempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -csh -pack=GeoPrimitives -version=GeoPrimitives-00-00-33-05 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/DetectorDescription  $* >${cmtGeoPrimitivestempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtGeoPrimitivestempfile}
  unset cmtGeoPrimitivestempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtGeoPrimitivestempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtGeoPrimitivestempfile}
unset cmtGeoPrimitivestempfile
exit $cmtcleanupstatus

