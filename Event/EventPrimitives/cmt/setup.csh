# echo "setup EventPrimitives EventPrimitives-00-00-45 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/cern.ch/sw/contrib/CMT/v1r25p20140131
endif
source ${CMTROOT}/mgr/setup.csh
set cmtEventPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if $status != 0 then
  set cmtEventPrimitivestempfile=/tmp/cmt.$$
endif
${CMTROOT}/${CMTBIN}/cmt.exe setup -csh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event  -no_cleanup $* >${cmtEventPrimitivestempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/${CMTBIN}/cmt.exe setup -csh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event  -no_cleanup $* >${cmtEventPrimitivestempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtEventPrimitivestempfile}
  unset cmtEventPrimitivestempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtEventPrimitivestempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtEventPrimitivestempfile}
unset cmtEventPrimitivestempfile
exit $cmtsetupstatus

