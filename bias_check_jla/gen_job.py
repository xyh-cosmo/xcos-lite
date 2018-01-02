import sys

if len(sys.argv) < 5:
    print 'usage:\n\tpython %s %s %s %s %s'%(sys.argv[0],'KEY_JOB_ID','KEY_LOGFILE','KEY_INI','JOB_FILENAME')
    sys.exit(0)

job_filename=sys.argv[4]
fp = open(job_filename,'w')

print >> fp, "#!/bin/bash"
print >> fp, "#$ -N %s"%(sys.argv[1])
print >> fp, "#$ -S /bin/bash"
print >> fp, "#$ -o %s"%(sys.argv[2])
print >> fp, "# join stderr and stdout\n"
  
print >> fp, "unset SGE_ROOT"
print >> fp, "cd /home/xyh/nfs_xyh/xcos-lite\n"

print >> fp, "#echo $NSLOTS"
print >> fp, "tmphosts=`cat $TMPDIR/machines`"
print >> fp, "hosts=`echo $tmphosts | sed 's/\ /,/g'`\n"


print >> fp, "mpirun -np $NSLOTS -H $hosts \\"
print >> fp, "-mca btl_openib_if_include mlx4_0:1 -mca btl sm,self,openib \\"
print >> fp, "-mca coll_fca_np 4 -mca coll_fca_enable 0 \\"
print >> fp, "bin/iSampler %s"%(sys.argv[3])
  
fp.close()


