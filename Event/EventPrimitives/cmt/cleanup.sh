# echo "cleanup EventPrimitives EventPrimitives-00-00-45 in /afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r25p20140131; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtEventPrimitivestempfile=`${CMTROOT}/${CMTBIN}/cmt.exe -quiet build temporary_name`
if test ! $? = 0 ; then cmtEventPrimitivestempfile=/tmp/cmt.$$; fi
${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event  $* >${cmtEventPrimitivestempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/${CMTBIN}/cmt.exe cleanup -sh -pack=EventPrimitives -version=EventPrimitives-00-00-45 -path=/afs/cern.ch/work/j/jhrdinka/FCC/ats-FCC/a-tracking-sw/Event  $* >${cmtEventPrimitivestempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtEventPrimitivestempfile}
  unset cmtEventPrimitivestempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtEventPrimitivestempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtEventPrimitivestempfile}
unset cmtEventPrimitivestempfile
return $cmtcleanupstatus

