Tips for Minimizing AWS costs
=============================

Switching instance types
------------------------
To run the inversion, you need to have an instance with significant compute power, but if you are only running the setup script or doing analysis of output data, then much of this compute power is going to waste. However, it is possible to switch the instance type of an existing instance if you expect to be doing less compute heavy tasks. See the `AWS Documentation on how to change your instance type <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-resize.html>`_.

Spot Instances
--------------
Spot instances take advantage of unused compute capacity on the AWS cloud, allowing users to launch instances at a 70-90% reduction in price compared to on-demand instances. However, this reduced pricing comes with the understanding that AWS can take back this extra capacity at any time, so your instance may be interrupted (with a 2 minute warning). Generally, interruptions are rare (~5% of instances get interrupted) and once interuppted your instance will be either Stopped or Terminated depending on your configuration. We recommend using spot instances for simulation that take hours (not days) as it can greatly reduce your EC2 costs. For information on how to launch a spot instance see `Create a Spot Instance Request <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/spot-requests.html#create-spot-instance-request-console-procedure>`_. For more information on how to avoid and handle interruptions checko out this post on `Best Practices <https://aws.amazon.com/blogs/compute/best-practices-for-handling-ec2-spot-instance-interruptions/>`_.

Allocating Storage Volume Size
------------------------------
You pay for the amount of storage space you allocate to an instance at creation time. The cost of storage is not high, but can become significant if you keep the storage volumes for long periods of time (weeks/months). It is best to only allocate the storage you need and delete the volume once finished with it to prevent these costs. Note that volumes can be dynamically expanded, but not easily shrunk, see the `AWS Documentation for details <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/requesting-ebs-volume-modifications.html>`_.

Exporting Data to S3
--------------------
Storing data in EBS volumes is more expensive than storing data in AWS' simple storage service (S3). Additionally, data costs in S3 are only charged on the amount of space you use, whereas EBS volumes charge you for the amount of space allocated.

For these reasons, it is best practice to, after running the IMI, push your needed output data to an S3 bucket for long term storage and usage, rather than storing the entire EBS volume. Resources for creating an s3 bucket and pushing data to it can be found here:

* `Creating an S3 Bucket <https://docs.aws.amazon.com/AmazonS3/latest/userguide/create-bucket-overview.html>`_
* `Uploading/Downloading files using the cp command <https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html#using-s3-commands-managing-objects-copy>`_
