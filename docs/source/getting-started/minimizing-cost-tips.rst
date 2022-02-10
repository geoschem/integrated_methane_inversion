Tips for minimizing AWS costs
=============================


Switching instance types
------------------------
The IMI can be very efficient when run on an EC2 instance with significant compute power. But if you only 
wish to run the IMI preview or analyze output data, then much of this compute power will be wasted. 

Thankfully, it is possible to switch the instance type of an existing instance if you expect to be doing less compute-heavy work. 
See the `AWS Documentation on how to change your instance type <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-resize.html>`_
for more information.


Spot instances
--------------
Spot instances take advantage of unused compute capacity on the AWS cloud, allowing users to launch instances at a 
70-90% reduction in price compared to on-demand instances. 

However, this reduced pricing comes with the understanding that AWS can take back this extra capacity at any time, 
so your instance may be interrupted (with a 2 minute warning), causing the IMI to crash. 
Interruptions are generally rare (~5% of instances get interrupted) and once interupted your instance will be either 
Stopped or Terminated depending on the EC2 configuration. 

We recommend using spot instances for inversions that take hours (not days) as it can greatly reduce your EC2 costs. 
For information on how to launch a spot instance see 
`Create a Spot Instance Request <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-requests.html#create-spot-instance-request-console-procedure>`_. 
For more information on how to avoid and handle interruptions check out this post on 
`Best Practices <https://aws.amazon.com/blogs/compute/best-practices-for-handling-ec2-spot-instance-interruptions/>`_.


.. _selectingStorageSize-label:

Selecting storage volume size
-----------------------------
AWS charges continuous fees for the storage volume provisioned to an EC2 instance. 
These fees can become significant if you retain the volume for long periods of time (weeks/months). 

It is best to provision only the amount of storage needed, and to delete your volume once finished with it to minimize costs. 

You can always add storage space after launching an EC2 instance, but it is very difficult to retroactively reduce storage space;
see the `AWS Documentation for details <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/requesting-ebs-volume-modifications.html>`_.

.. note::
  When unsure of the storage needs for an inversion, we recommend starting small. A good starting point is ~100 GB. 
  
  To determine your true storage needs, first ``ssh`` into the instance and run a 1-week inversion for 
  your region of interest. When the 1-week inversion is complete, check how much storage has been used. 
  From there, you can scale-up the storage according to your actual period of interest. Consider that
  the AMI itself takes about 20 GB of storage.

  For example, if after the 1-week inversion you find that 75/100 GB are occupied, then you should budget
  75 - 20 = 55 GB per inversion week. If you want to perform a 1-year inversion, then increasing the storage 
  to 3.5 TB will leave you with about 500 GB of additional space to work with once the inversion is complete.


.. _exportingS3-label:

Exporting data to S3
--------------------
Storing data in EBS volumes is more expensive than storing data in Amazon S3. 
Additionally, with S3 you are only charged for the amount of space you use, whereas EBS volumes 
charge you for the amount of space provisioned.

For these reasons, after running the IMI, we recommend pushing your output data to an S3 bucket 
for long term storage, rather than retaining the entire EBS volume. 
Resources for creating an S3 bucket and pushing data to it can be found here:

* `Creating an S3 Bucket <https://docs.aws.amazon.com/AmazonS3/latest/userguide/create-bucket-overview.html>`_
* `Uploading/Downloading files using the cp command <https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html#using-s3-commands-managing-objects-copy>`_