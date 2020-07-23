
## Stuff to get started
'''
https://www.tensorflow.org/tutorials/quickstart/beginner
'''
import os
#os.environ["CUDA_VISIBLE_DEVICES"]="-1"
import tensorflow as tf

# Check available devices

visible_devices = tf.config.get_visible_devices(
    device_type=None
)
print(visible_devices)

# Tell tensorflow, that we run on CPU
#sess = tf.Session(config=tf.ConfigProto(device_count={'GPU': 0}))

mnist = tf.keras.datasets.mnist

(x_train, y_train), (x_test, y_test) = mnist.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0

model = tf.keras.models.Sequential([
  tf.keras.layers.Flatten(input_shape=(28, 28)),
  tf.keras.layers.Dense(128, activation='relu'),
  tf.keras.layers.Dropout(0.2),
  tf.keras.layers.Dense(10)
])

# model execution
predictions = model(x_train[:1]).numpy()
print(predictions)

print(tf.nn.softmax(predictions).numpy())

# define loss
loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

print(loss_fn(y_train[:1], predictions).numpy())


# compile model
model.compile(optimizer='adam',
              loss=loss_fn,
              metrics=['accuracy'])

# Train
model.fit(x_train, y_train, epochs=5)

model.evaluate(x_test,  y_test, verbose=2)


probability_model = tf.keras.Sequential([
  model,
  tf.keras.layers.Softmax()
])

probability_model(x_test[:5])

